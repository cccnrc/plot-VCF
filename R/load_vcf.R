#' Load a VCF file
#'
#' This function loads a VCF file as a matrix.
#'
#' @param VCF_FILE Path to the input VCF file
#' @param ASSEMBLY (optional) assembly of the VCF file (default hg38)
#' @param SAMPLE (optional) samples to plot
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants, default is just position
#' @return A matrix of the infile
load_vcf <- function( VCF_FILE, ASSEMBLY="hg38", SAMPLE="ALL", VAR_FLAG="POS" ){

  VCF_DATA <- suppressWarnings(VariantAnnotation::readVcf( VCF_FILE, ASSEMBLY ))
  POS_COL <- data.frame(rowRanges(VCF_DATA)[,"paramRangeID"])[,c('seqnames', 'start')]
  VAR_COL <- data.frame(rowRanges(VCF_DATA))[,c('QUAL', 'FILTER')]

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
    ### extract AF column if present
    if ( 'AF' %in% names(geno(VCF_DATA)) ) {
      AF_COL <- geno(VCF_DATA)$AF[,SAMPLE]
    }
    ### extract DP column if present
    if ( 'DP' %in% names(geno(VCF_DATA)) ) {
      DP_COL <- geno(VCF_DATA)$DP[,SAMPLE]
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
    ### extract AF column if present
    if ( 'AF' %in% names(geno(VCF_DATA)) ) {
      AF_COL <- geno(VCF_DATA)$AF
    }
    ### extract DP column if present
    if ( 'DP' %in% names(geno(VCF_DATA)) ) {
      DP_COL <- geno(VCF_DATA)$DP
    }
  }

  ### if AF column found
  if ( exists("AF_COL") ) {
    if (( length(SAMPLE) > 1 )||( SAMPLE == "ALL" )) {
      AF_MED <- apply(AF_COL,1,function(v) median(as.numeric(v),na.rm = T))
      AF_MEDF <- data.frame( 'AF' = AF_MED )
    } else {
      AF_MEDF <- data.frame( 'AF' = as.numeric(AF_COL) )
    }
  }

  ### if DP column found
  if ( exists("DP_COL") ) {
    if (( length(SAMPLE) > 1 )||( SAMPLE == "ALL" )) {
      DP_MED <- apply(DP_COL,1,function(v) median(as.numeric(v),na.rm = T))
      DP_MEDF <- data.frame( 'DP' = DP_MED )
    } else {
      DP_MEDF <- data.frame( 'DP' = as.numeric(DP_COL) )
    }
  }


  ### if AF found add it to the dataframe
  if ( exists("AF_MEDF") ) {
    VCF_BODY_PRE <- cbind( POS_COL, rownames(GENO_COL), VAR_COL, AF_MEDF )
  } else {
    VCF_BODY_PRE <- cbind( POS_COL, rownames(GENO_COL), VAR_COL )
  }

  ### if DP found add it to the dataframe
  if ( exists("DP_MEDF") ) {
    VCF_BODY_PRE <- cbind( VCF_BODY_PRE, DP_MEDF )
  }

  ### add GT columns
  VCF_BODY <- cbind( VCF_BODY_PRE, GENO_COL )

  ### put correct colnames
  colnames( VCF_BODY ) <- c( c( 'CHROM', 'POS', 'VAR' ), colnames(VAR_COL), colnames(VCF_BODY)[6:length(colnames(VCF_BODY))] )

  ### check VAR_FLAG was found (if specified)
  if ( VAR_FLAG != "POS" ) {
    if ( VAR_FLAG %in% colnames(VCF_BODY) ) {
      cat( "    -> VAR_FLAG column found and loaded\n" )
    } else {
      cat( "\n\t-> ERROR: cannot find VAR_FLAG value in VCF column!\n" )
      stop()
    }
  }

  VCF_BODY
}
