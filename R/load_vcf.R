#' Load a VCF file
#'
#' This function loads a VCF file as a matrix.
#'
#' @param VCF_FILE Path to the input VCF file
#' @param SEQINFO the object with chromosome characteristics: load_chr()
#' @param ASSEMBLY (optional) assembly of the VCF file (default hg38)
#' @param SAMPLE (optional) samples to plot (default is all samples in VCF file)
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants, default is just position
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return A matrix of the infile
load_vcf <- function( VCF_FILE, SEQINFO, ASSEMBLY="hg38", SAMPLE="ALL", VAR_FLAG="POS", CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY") ){
  ### prepare VCF_FILE for reading
  getExtension <- function(FILE_PATH) {
      ex <- strsplit(basename(FILE_PATH), split="\\.")[[1]]
      cat('  -> passed VCF file extension:', ex[length(ex)], '\n')
      return(ex[length(ex)])
  }
  compressVCF <- function( VCF_PATH ) {
    cat('  -> compressing VCF file ...\n')
    VCF_GZ <- Rsamtools::bgzip( VCF_PATH, tempfile() )
    return( VCF_GZ )
  }
  indexVCFGZ <- function( VCF_GZ ) {
    cat('  -> indexing VCF file ...\n')
    idx <- Rsamtools::indexTabix( VCF_GZ, "vcf" )
    tab <- Rsamtools::TabixFile( VCF_GZ, idx )
    return( tab )
  }
  getVCFchr <- function ( VCF_GZ ) {
    CHR_VCF_GZ <- Rsamtools::headerTabix( VCF_GZ, "vcf" )$seqnames
    return(CHR_VCF_GZ)
  }
  getPARAM <- function( CHR_NAMES, SEQINFO, SAMPLE ) {
    SEQ_CHR <- SEQINFO[ CHR_NAMES ]
    GR_CHR <- GenomicRanges::GRanges(SEQ_CHR)
    PARAM <- VariantAnnotation::ScanVcfParam( which = GR_CHR )
    if (( length(SAMPLE) > 1 ) || ( SAMPLE != "ALL" )) {
      PARAM <- VariantAnnotation::ScanVcfParam( which = GR_CHR, samples = SAMPLE )
    }
    return( PARAM )
  }
  calculateAF <- function( VCF_DATA, SAMPLE ) {
    GH <- data.frame(geno(header(VCF_DATA)))
    GHE <- GH[ GH[,'Type'] == "Integer" & GH[,'Number'] == 1, ]
    GHE_FLOAT <- GH[ GH[,'Type'] == "Float", ]
    if ( ( "AD" %in% rownames(GH) ) && ( "DP" %in% rownames(GHE) ) ) {
      if (( length(SAMPLE) > 1 ) || ( SAMPLE != "ALL" )) {
        if ( length(SAMPLE) > 1 ) {
          TMP_DB_AD <- data.frame(geno(VCF_DATA)[["AD"]])[,SAMPLE]
          TMP_DB_DP <- data.frame(geno(VCF_DATA)[["DP"]])[,SAMPLE]
        } else {
          AD_VEC <- vector()
          for ( VAL in data.frame(geno(VCF_DATA)[["AD"]])[,SAMPLE] ) {
            AD_VEC <- c(AD_VEC, VAL[length(VAL)])
          }
          TMP_DB_AD <- data.frame( "SAMPLE" = AD_VEC )
          TMP_DB_DP <- data.frame( "DP" = data.frame(geno(VCF_DATA)[["DP"]])[,SAMPLE] )
        }
      } else {
        TMP_DB_AD <- data.frame(geno(VCF_DATA)[["AD"]])
        TMP_DB_DP <- data.frame(geno(VCF_DATA)[["DP"]])
      }
      AF_LIST <- list()
      for ( c in 1:ncol( TMP_DB_AD ) )
      {
        SAMPLE_AF <- vector()
        for ( i in 1:nrow(TMP_DB_AD) )
        {
          AD_VALUES <- TMP_DB_AD[i,c][[1]]
          AD_ALT <- as.numeric(AD_VALUES[ length(AD_VALUES) ])
          AF_VALUE <- AD_ALT/as.numeric(TMP_DB_DP[i,c])
          SAMPLE_AF <- c(SAMPLE_AF, AF_VALUE)
        }
        AF_LIST[[c]] <- SAMPLE_AF
      }
      names( AF_LIST ) <- colnames( TMP_DB_AD )
      TMP_DB_AF <- data.frame(AF_LIST)
      if (( length(SAMPLE) > 1 ) || ( SAMPLE == "ALL" )) {
        TMP_COL_AF <- apply(TMP_DB_AF,1,function(v) median(as.numeric(v),na.rm = T))
      } else {
        TMP_COL_AF <- as.numeric(TMP_DB_AF[,1])
      }
      VAR_COL <- data.frame( "AF" = as.numeric(TMP_COL_AF) )
      cat( "      -> VAR_FLAG column: AF calculated from AD/DP ...\n" )
    } else {
      stop("      -> can not calculate AF! AD and/or DP absent from VCF file ...")
    }
    return(VAR_COL)
  }

  ### compress and index file (if not done yet)
  VCF_EXT <- getExtension(VCF_FILE)
  if (( VCF_EXT == "vcf" ) || ( VCF_EXT == "VCF" )) {
    #cat(' -> passed file is .vcf Compressing and indexing it ...\n')
    VCF_GZ <- compressVCF( VCF_FILE )
    VCF_TAB <- indexVCFGZ( VCF_GZ )
    ### check and get chromosmes actually in VCF file
    VCF_CHR_NAMES <- getVCFchr( VCF_GZ )
  } else if (( VCF_EXT == "gz" ) || ( VCF_EXT == "GZ" )) {
    #cat(' -> passed file is .vcf.gz Indexing it ...\n')
    VCF_TAB <- indexVCFGZ( VCF_FILE )
    ### check and get chromosmes actually in VCF file
    VCF_CHR_NAMES <- getVCFchr( VCF_FILE )
  } else {
    cat('\n')
    cat(' -> passed VCF file is neither .vcf or .vcf.gz! (see documentation)...\n')
    stop()
  }
  ### compare CHR_NAMES with chromosomes actually in VCF file
  CHR_ABS <- setdiff( CHR_NAMES, VCF_CHR_NAMES )
  if ( length( CHR_ABS ) > 0 ) {
    for ( CHR_AB in CHR_ABS ) {
      cat( "    -> VCF, no variant found in:", CHR_AB,"(will remove) \n" )
    }
    CHR_NAMES <- CHR_NAMES[!CHR_NAMES %in%CHR_ABS]
  }
  ### adapt SEQINFO
  SEQINFO <- SEQINFO[ CHR_NAMES ]
  ### actually read and load VCF file
  VCF_PARAM <- getPARAM( CHR_NAMES, SEQINFO, SAMPLE )
  cat('  -> reading VCF file ...\n')
  VCF_DATA <- suppressWarnings(VariantAnnotation::readVcf( VCF_TAB, ASSEMBLY, VCF_PARAM ))
  POS_COL <- data.frame(rowRanges(VCF_DATA)[,"paramRangeID"])[,c('seqnames', 'start')]
  VCF_SAMPLES <- rownames(colData( VCF_DATA ))
  cat( "    -> VCF samples:", length(VCF_SAMPLES),"\n" )
  ### check all specified samples are in VCF
  if (( length(SAMPLE) > 1 ) || ( SAMPLE != "ALL" )) {
    for ( SAM in SAMPLE )
    {
      if ( ! SAM %in% VCF_SAMPLES ) {
        cat( "\t-> ERROR: cannot find specified sample:", SAM, "in VCF!\n" )
        stop()
      }
    }
    cat( "    -> all specified samples", paste("(n. ", length(SAMPLE), ")", sep = ''), "in VCF!\n" )
  }
  ### if user wants POS do not need to store the rest of database
  if ( VAR_FLAG != "POS" ) {
    if ( VAR_FLAG == "QUAL" ) {
      TMP_COL_DB <- data.frame(rowRanges(VCF_DATA))[,'QUAL']
      VAR_COL <- data.frame( "QUAL" = TMP_COL_DB )
    } else {
      ### select genotype columns that are numeric with single data
      GH <- data.frame(geno(header(VCF_DATA)))
      GHE <- GH[ GH[,'Type'] == "Integer" & GH[,'Number'] == 1, ]
      GHE_FLOAT <- GH[ GH[,'Type'] == "Float", ]
      if (( VAR_FLAG %in% rownames(GHE) )||( VAR_FLAG %in% rownames(GHE_FLOAT) )) {
        ### get flag number
        NAME <- VAR_FLAG
        if (( length(SAMPLE) > 1 ) || ( SAMPLE != "ALL" )) {
          TMP_DB <- data.frame(geno(VCF_DATA)[[NAME]])[,SAMPLE]
        } else {
          TMP_DB <- data.frame(geno(VCF_DATA)[[NAME]])
        }
        if (( length(SAMPLE) > 1 ) || ( SAMPLE == "ALL" )) {
          TMP_COL <- apply(TMP_DB,1,function(v) median(as.numeric(v),na.rm = T))
        } else {
          TMP_COL <- as.numeric(TMP_DB)
        }
        TMP_COL_DB <- data.frame( NAME = TMP_COL )
        colnames(TMP_COL_DB) <- NAME
        VAR_COL <- TMP_COL_DB
        cat( "      -> VAR_FLAG column:", NAME, "found and loaded...\n" )
      } else {
        ### if VAR_FLAG is AF and it is absent from VCF calculate it from AD adn DP
        if ( VAR_FLAG == "AF" ) {
          VAR_COL <- calculateAF( VCF_DATA, SAMPLE )
        }
      }
    }
  }

  ### get GT value
  if (( length(SAMPLE) > 1 ) || ( SAMPLE != "ALL" )) {
    ### extract GT column
    if ( 'GT' %in% names(geno(VCF_DATA)) ) {
      GENO_COL <- geno(VCF_DATA)$GT[,SAMPLE]
      if ( length(SAMPLE) == 1 ) {
        GENO_COL <- data.frame(SAMPLE = GENO_COL)
        colnames(GENO_COL) <- SAMPLE
      }
    } else {
      cat( "\n\t-> ERROR: cannot find GT value in VCF genotype column!\n" )
      stop()
    }
  ### if SAMPLE not specified load all samples
  } else {
    ### extract GT column
    if ( 'GT' %in% names(geno(VCF_DATA)) ) {
      GENO_COL <- geno(VCF_DATA)$GT
    } else {
      cat( "\n\t-> ERROR: cannot find GT value in VCF genotype column!\n" )
      stop()
    }
  }

  ### add GT columns
  VCF_BODY <- cbind( POS_COL, rownames(GENO_COL), VAR_COL, GENO_COL )

  ### put correct colnames
  colnames( VCF_BODY ) <- c( c( 'CHROM', 'POS', 'VAR' ), colnames(VAR_COL), colnames(GENO_COL) )

  ### check VAR_FLAG was found (if specified)
  if ( VAR_FLAG != "POS" ) {
    if ( ! VAR_FLAG %in% colnames(VCF_BODY) ) {
      cat( "\n\t-> ERROR: cannot find VAR_FLAG value in VCF column!\n" )
      stop()
    }
  }
  ### return correct objects
  RETURN_LIST <- list("VCF" = VCF_BODY, "SEQINFO" = SEQINFO)
  return( RETURN_LIST )
}
