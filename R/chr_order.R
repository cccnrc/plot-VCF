#' return a vector of chromosome names in proper genomic order
#'
#' This function takes in a vector of chromosome names and return them in proper order
#'
#' @param CHR_VECTOR vector with chromosme names
#' @return ordered vector of CHR_VECTOR
chr_order <- function( CHR_VECTOR ){
  CHR_NAMES38 <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  CHR_NAMES37 <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  ### check which assembly is the passed data
  CHR_VECTOR <- unique(CHR_VECTOR)
  OUTPUT_VECTOR <- vector()
  if ( all(CHR_VECTOR %in% CHR_NAMES38) ) {
    for (CHR in CHR_NAMES38)
    {
      if ( CHR %in% CHR_VECTOR ) {
        OUTPUT_VECTOR <- c(OUTPUT_VECTOR, CHR)
      }
    }
  } else if ( all(CHR_VECTOR %in% CHR_NAMES37) ) {
    for (CHR in CHR_NAMES37)
    {
      if ( CHR %in% CHR_VECTOR ) {
        OUTPUT_VECTOR <- c(OUTPUT_VECTOR, CHR)
      }
    }
  } else {
    stop("    -> unrecognized passed chromosome names:", CHR_VECTOR )
  }
  return( OUTPUT_VECTOR )
}
