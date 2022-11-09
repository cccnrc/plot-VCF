#' Adapt samples specified in COLOR_SAMPLE
#'
#' This function iterates through COLOR_SAMPLE and checks it
#' It the returns a little DF of 2 columns with each specified sample and corresponding group
#'
#' @param COLOR_SAMPLE list with groups and names of sample to color in the final plot
#' @param SAMPLE samples to plot (default is all samples in VCF file)
#' @return A 2-columns dataframe with sample-group correspondence
adapt_color_sample <- function( COLOR_SAMPLE, SAMPLE ){
  ### check COLOR_SAMPLE if passed
  COLOR_SAMPLE_GROUPS <- length( COLOR_SAMPLE )
  if ( COLOR_SAMPLE_GROUPS > 1 ) {
    ### check they do not specify both SAMPLE and COLOR_SAMPLE
    if ( ( length(SAMPLE) > 1 ) || ( SAMPLE != "ALL" ) ) {
      warning(' -> both SAMPLE and COLOR_SAMPLE specified. Only COLOR_SAMPLE will be used ...')
    }
    # cat("    ->", COLOR_SAMPLE_GROUPS, "groups passed to COLOR_SAMPLE\n")
    COLOR_SAMPLE_VEC <- vector()
    COLOR_SAMPLE_GROUP_VEC <- vector()
    for ( i in 1:COLOR_SAMPLE_GROUPS ) {
      # cat("      -> group:", names(COLOR_SAMPLE[i]), "- n:", length(COLOR_SAMPLE[[i]]), "\n")
      ### stop if no samples in this COLOR_SAMPLE group
      if ( length(COLOR_SAMPLE[[i]]) == 0 ) {
        stop(" => no samples specified in", COLOR_SAMPLE_GROUPS[i], "COLOR_SAMPLE group")
      }
      ### put all values in a vector
      COLOR_SAMPLE_VEC <- c( COLOR_SAMPLE_VEC, COLOR_SAMPLE[[i]] )
      COLOR_SAMPLE_GROUP_VEC <- c( COLOR_SAMPLE_GROUP_VEC, names( COLOR_SAMPLE[i]) )
    }
    ### remove duplicates in COLOR_SAMPLE
    COLOR_SAMPLE_DUP <- COLOR_SAMPLE_VEC[duplicated(COLOR_SAMPLE_VEC)]
    if ( length(COLOR_SAMPLE_DUP) > 0 ) {
      stop(" => duplicated samples passed in COLOR_SAMPLE!")
    }
    ### if specified, restrict analysis on samples in COLOR_SAMPLE
    NEW_SAMPLE <- COLOR_SAMPLE_VEC
  }
  COLOR_SAMPLE_DF <- data.frame( "SAM" = NEW_SAMPLE, "GROUP" = COLOR_SAMPLE_GROUP_VEC )
  return( COLOR_SAMPLE_DF )
}
