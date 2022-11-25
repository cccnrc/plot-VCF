#' return a vector with null GT removed
#'
#' This function takes in a GT vector object and returns it without null GTs
#'
#' @param GT_VECTOR vector object with GT (as extracted by VariantAnnotation::readVcf())
#' @param GT_NULL (optional) vector passing GT values to consider null (default: c( '0/0', './.', '0|0', '.|.', '0', '.' ))
#' @return non null GT vector
gt_null <- function( GT_VECTOR, GT_NULL=c( '0/0', './.', '0|0', '.|.', '0', '.' ) ){
  GT_N <- vector()
  for ( VALUE in GT_VECTOR )
  {
    if ( ! VALUE %in% GT_NULL ) {
      GT_N <- c( GT_N, VALUE )
    }
  }
  return( GT_N )
}
