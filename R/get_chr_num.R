#' Extract number of variants for each chromosome in the VCF file
#'
#' This function extract the number of variants for each chromosome from the modeled VCF file
#'
#' @param SAMPLE_GT_CLEAN the modeled VCF file
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return A vector with number of variants in the modeled VCF file for each specified chrosomome
get_chr_num <- function(SAMPLE_GT_CLEAN, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")){
  CHR_N <- vector()
  for ( CHR in CHR_NAMES )
  {
    N <- nrow( SAMPLE_GT_CLEAN[ SAMPLE_GT_CLEAN$CHROM == CHR, ] )
    # cat('    - ', CHR, ': ', N, '\n'  )
    CHR_N <- c( CHR_N, N )
  }
  CHR_N
}
