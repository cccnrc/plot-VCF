#' Create a vector with height to plot each variant
#'
#' This function defines the Y axis coordinate to plot each variant, based on the total variants in the VCF file
#'
#' @param CHR_N the vector with number of variants for each chrosomome (as returned by get_chr_num() function)
#' @param ORDERED (optional) if user want to have ordered variant height in the plot or random
#' @return A vector with height to plot each variant
get_chr_plot_step <- function(CHR_N, ORDERED=FALSE){
  if (( ORDERED != FALSE )&( ORDERED != F )) {
    CHR_N_MAX <- max( CHR_N )
    STEP <- 1/CHR_N_MAX
    VAR_HEIGHT <- vector()
    for ( TOT in CHR_N )
    {
      if ( TOT > 0 ) {
        for ( i in 1:TOT )
        {
          HEIGHT_I <- STEP * i
          VAR_HEIGHT <- c( VAR_HEIGHT, HEIGHT_I )
        }
      }
    }
    CHR_Y <- VAR_HEIGHT
  }
  else {
    ### if un-ordered is requested it gives random Y axis values
    GL <- runif(n=sum(CHR_N), min=1e-12, max=.99)
    CHR_Y <- GL
  }
  return(CHR_Y)
}
