#' Load a VCF file
#'
#' This function loads a VCF file as a matrix.
#'
#' @param VCF_FILE Path to the input VCF file
#' @param ASSEMBLY (optional) assembly of the VCF file (default hg38)
#' @param SAMPLE (optional) samples to plot (default is all samples in VCF file)
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants, default is just position
#' @return A matrix of the infile
load_vcf <- function( VCF_FILE, ASSEMBLY="hg38", SAMPLE="ALL", VAR_FLAG="POS" ){

  VCF_DATA <- suppressWarnings(VariantAnnotation::readVcf( VCF_FILE, ASSEMBLY ))
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
      if ( VAR_FLAG %in% rownames(GHE) ) {
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

  VCF_BODY
}
