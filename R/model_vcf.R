#' Modify a VCF file to get a row for each variant for each sample and a GT column
#'
#' This function modify a loaded VCF (load_vcf) file.
#'
#' @param VCF_DATA the loaded VCF file
#' @return modified dataframe
model_vcf <- function(VCF_DATA){
  VCF_COLNAMES <- colnames(VCF_DATA)
  ### extract sample names from VCF_DATA
  SAMPLE_NAMES <- VCF_COLNAMES[5:length(VCF_COLNAMES)]
  GENO_DATA <- VCF_DATA[, 1:4]
  ### create the first sample column
  SAMPLE <- SAMPLE_NAMES[1]
  SAMPLE_COL <- rep( SAMPLE, nrow( GENO_DATA ) )
  SAMPLE_DATA <- cbind( SAMPLE_COL, GENO_DATA, VCF_DATA[ ,SAMPLE ] )
  SAMPLE_DATA_COLNAMES <- colnames(SAMPLE_DATA)
  SAMPLE_DATA_COLNAMES[1] <- 'IND'
  SAMPLE_DATA_COLNAMES[2] <- 'CHROM'
  SAMPLE_DATA_COLNAMES[length(SAMPLE_DATA_COLNAMES)] <- 'GT'
  colnames(SAMPLE_DATA) <- SAMPLE_DATA_COLNAMES
  SAMPLE_GT <- SAMPLE_DATA
  ### add the rest of samples (if present)
  if ( length(SAMPLE_NAMES) > 1 ) {
    for (SAMPLE in SAMPLE_NAMES[2:length(SAMPLE_NAMES)])
    {
      SAMPLE_COL <- rep( SAMPLE, nrow( GENO_DATA ) )
      SAMPLE_DATA <- cbind( SAMPLE_COL, GENO_DATA, VCF_DATA[ ,SAMPLE ] )
      colnames(SAMPLE_DATA) <- SAMPLE_DATA_COLNAMES
      SAMPLE_GT <- rbind( SAMPLE_GT, SAMPLE_DATA )
    }
  }

  ### remove rows with null GT
  NULL_GT <- c( './.', '0/0', '.|.', '0|0' )
  SAMPLE_GT_CLEAN <- SAMPLE_GT[ SAMPLE_GT[, 'GT'] != '.', ]
  for (NG in NULL_GT)
  {
    SAMPLE_GT_CLEAN <- SAMPLE_GT_CLEAN[!startsWith( SAMPLE_GT_CLEAN[,'GT'], NG ),]
  }
  SAMPLE_GT_CLEAN
}
