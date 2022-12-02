#' Load a VCF file
#'
#' This function loads a VCF file as a matrix.
#'
#' @param VCF_DATA loaded VCF file (read_vcf() output)
#' @param SEQINFO the object with chromosome characteristics: load_chr()
#' @param SAMPLE (optional) samples to plot (default is all samples in VCF file)
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants, default is just position
#' @return A matrix of the infile
load_vcf <- function( VCF_DATA, SEQINFO, SAMPLE="ALL", VAR_FLAG="POS" ){
  ### derive AF value from AD and DP, if AF flag absent from VCF file
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

  ### model the VCF file to export needed fields
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
  } else {
    VAR_COL <- data.frame( "start" = POS_COL$start )
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
