#' Create the plot based on modeled VCF
#'
#' This function defines the Y axis coordinate to plot each variant, based on the total variants in the VCF file
#'
#' @param MODELED_VCF the database from the VCF file (load_vcf() and model_vcf())
#' @param SEQINFO the object with chromosome characteristics: load_chr()
#' @param VAR_Y (optional) the Y-axis coordinates for the variants
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants
#' @param SPACELINE (optional) show line between chromosomes
#' @param THRESHOLD (optional) a threshold line to use as Y-axis
#' @param XLIM (optional) limits for Y-axis
#' @param SHAPE (optional) different shape for each sample
#' @param SAMPLE_DB (optional) 2-column dataframe with sample-group to color plot accordingly
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return The VCF variant plot
make_plot <- function(MODELED_VCF, SEQINFO, VAR_Y=FALSE, VAR_FLAG="POS", SPACELINE=FALSE, SHAPE=FALSE, SAMPLE_DB=FALSE, THRESHOLD=FALSE, XLIM=FALSE, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")){
  ### set Y-axis limits as specified by user
  if ( XLIM != FALSE ) {
    if ( length(XLIM) == 1 ) {
      XLIMITS <- c(XLIM, NA)
    } else {
      XLIMITS <- XLIM
    }
  } else {
    XLIMITS <- c(NA, NA)
  }
  ### if user specified a column to use to plot data use this
  Y_AXIS_TEXT <- ggplot2::element_blank()
  if ( VAR_FLAG != "POS" ) {
    if ( VAR_FLAG %in% colnames(MODELED_VCF) ){
      MODELED_VCF  <- MODELED_VCF[ MODELED_VCF[,'CHROM'] %in% CHR_NAMES, c('IND','CHROM','POS', VAR_FLAG)]
      VAR_Y <- MODELED_VCF[,VAR_FLAG]
      Y_AXIS_TEXT <- ggplot2::element_text()
    }
    else {
      cat( "\n\t-> ERROR: undefined VCF flag selected with VAR_FLAG!\n" )
      stop()
    }
  }
  else {
    MODELED_VCF  <- MODELED_VCF[ MODELED_VCF[,'CHROM'] %in% CHR_NAMES, c('IND','CHROM','POS')]
  }
  ### default if no THRESHOLD was specified
  Y_THRESHOLD <- 0
  Y_THRESHOLD_THICKNESS <- 0
  if ( THRESHOLD != FALSE ) {
    Y_THRESHOLD <- THRESHOLD
    Y_THRESHOLD_THICKNESS <- 0.4
  }
  ### if specified, create factors to plot samples in different colors
  if ( (exists("SAMPLE_DB")) && ( is.data.frame(get('SAMPLE_DB')) ) ) {
    COLOR_SAMPLE_VEC <- vector()
    for ( i in 1:nrow(MODELED_VCF) )
    {
      COLOR_SAMPLE_VEC <- c( COLOR_SAMPLE_VEC, SAMPLE_DB[ SAMPLE_DB$SAM == MODELED_VCF$IND[i], 'GROUP' ] )
    }
    cat( '    -> adding COLOR_SAMPLE column to VCF plot ... \n' )
    COLOR_SAMPLE_DF <- data.frame( 'group' = COLOR_SAMPLE_VEC )
    MODELED_VCF <- cbind( MODELED_VCF, COLOR_SAMPLE_DF )
  }
  ### change plot color accordingly
  if ( "group" %in% colnames(MODELED_VCF) ) {
    MODELED_VCF$group <- factor( MODELED_VCF$group )
    PLOT_COLOR <- "group"
    LEGEND_COLOR <- "legend"
  } else {
    PLOT_COLOR <- "COL"
    LEGEND_COLOR <- "none"
  }
  ### order the DF based on chromosomes
  MODELED_VCF <- MODELED_VCF %>% arrange(factor( CHROM, levels = CHR_NAMES ))
  gr_geno <- suppressWarnings(GenomicRanges::makeGRangesFromDataFrame(MODELED_VCF, keep.extra.columns = TRUE,
                                      ignore.strand = TRUE, seqinfo = SEQINFO,
                                      seqnames.field = "CHROM", start.field = "POS",
                                      end.field = "POS"))
  ### specify color palette
  COL <- viridis::viridis(length(CHR_NAMES))
  ### specify shape palette
  gr_geno$IND <- factor( gr_geno$IND )
  if ( SHAPE != FALSE ) {
    if ( nlevels(gr_geno$IND) < 6 ) {
      SHAPE_SCALE <- 15:(16+nlevels(gr_geno$IND))
    } else {
      SHAPE_SCALE <- 1:nlevels(gr_geno$IND)
    }
    LEGEND_SHAPE <- "legend"
  } else {
    SHAPE_SCALE <- rep(19, nlevels(gr_geno$IND))
    LEGEND_SHAPE <- "none"
  }
  ### actually create the plot
  if ( "group" %in% colnames(MODELED_VCF) ) {
    VCF_PLOT <- suppressWarnings(suppressMessages(ggbio::plotGrandLinear(gr_geno, ggplot2::aes( y = VAR_Y, shape = IND, color = group  ),
                        space.skip = 0.01,
                        xlab = "Chromosome",
                        ylab = VAR_FLAG,
                        spaceline = SPACELINE ) +
                        ggplot2::theme_bw() +
                        ggplot2::scale_shape_manual(values=SHAPE_SCALE, name = '') +
                        ggplot2::scale_y_continuous(limits = XLIMITS) +
                        ggplot2::guides( colour = LEGEND_COLOR, shape = LEGEND_SHAPE ) +
                        ggplot2::geom_hline( yintercept=Y_THRESHOLD, linetype="longdash", color = "red", lwd = Y_THRESHOLD_THICKNESS) +
                        ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                              panel.background=ggplot2::element_blank(),
                              panel.border=ggplot2::element_blank(),
                              axis.line=ggplot2::element_blank(),
                              axis.text.y=Y_AXIS_TEXT,
                              axis.title.x=ggplot2::element_blank(),
                              axis.title.y=Y_AXIS_TEXT,
                              legend.position = "top",
                              legend.title=ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank(),
                              panel.grid.major=ggplot2::element_blank())))
  } else {
    VCF_PLOT <- suppressWarnings(suppressMessages(ggbio::plotGrandLinear(gr_geno, ggplot2::aes( y = VAR_Y, shape = IND  ),
                        space.skip = 0.01,
                        xlab = "Chromosome",
                        ylab = VAR_FLAG,
                        color = COL,
                        spaceline = SPACELINE ) +
                        ggplot2::theme_bw() +
                        ggplot2::scale_shape_manual(values=SHAPE_SCALE, name = '') +
                        ggplot2::scale_y_continuous(limits = XLIMITS) +
                        ggplot2::guides( colour = LEGEND_COLOR, shape = LEGEND_SHAPE ) +
                        ggplot2::geom_hline( yintercept=Y_THRESHOLD, linetype="longdash", color = "red", lwd = Y_THRESHOLD_THICKNESS) +
                        ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                              panel.background=ggplot2::element_blank(),
                              panel.border=ggplot2::element_blank(),
                              axis.line=ggplot2::element_blank(),
                              axis.text.y=Y_AXIS_TEXT,
                              axis.title.x=ggplot2::element_blank(),
                              axis.title.y=Y_AXIS_TEXT,
                              legend.position = "top",
                              panel.grid.minor.y = ggplot2::element_blank(),
                              panel.grid.major=ggplot2::element_blank())))
  }
  cat('\n')
  VCF_PLOT
}
