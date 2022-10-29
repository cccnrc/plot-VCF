#' Create the plot based on modeled VCF
#'
#' This function defines the Y axis coordinate to plot each variant, based on the total variants in the VCF file
#'
#' @param MODELED_VCF the database from the VCF file (load_vcf() and model_vcf())
#' @param SEQINFO the object with chromosome characteristics: load_chr()
#' @param VAR_Y (optional) the Y-axis coordinates for the variants
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants
#' @param THRESHOLD (optional) a threshold line to use as Y-axis
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return The VCF variant plot
make_plot <- function(MODELED_VCF, SEQINFO, VAR_Y=FALSE, VAR_FLAG="POS", THRESHOLD=FALSE, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")){
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
  ### order the DF based on chromosomes
  MODELED_VCF <- MODELED_VCF %>% arrange(factor( CHROM, levels = CHR_NAMES ))
  gr_geno <- suppressWarnings(GenomicRanges::makeGRangesFromDataFrame(MODELED_VCF, keep.extra.columns = TRUE,
                                      ignore.strand = TRUE, seqinfo = SEQINFO,
                                      seqnames.field = "CHROM", start.field = "POS",
                                      end.field = "POS"))
  COL <- viridis::viridis(length(CHR_NAMES))
  VCF_PLOT <- suppressWarnings(suppressMessages(ggbio::plotGrandLinear(gr_geno, ggplot2::aes(y = VAR_Y),
                      space.skip = 0.01,
                      xlab = "Chromosome",
                      ylab = "",
                      color = COL) +
                      ggplot2::theme_bw() +
                      # scale_x_discrete(labels = seqinfo) +
                      ggplot2::geom_hline( yintercept=Y_THRESHOLD, linetype="longdash", color = "red", lwd = Y_THRESHOLD_THICKNESS) +
                      ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                            panel.background=ggplot2::element_blank(),
                            panel.border=ggplot2::element_blank(),
                            axis.line=ggplot2::element_blank(),
                            axis.text.y=Y_AXIS_TEXT,
                            axis.title.x=ggplot2::element_blank(),
                            axis.title.y=Y_AXIS_TEXT,
                            legend.position = "none",
                            #panel.grid.minor=ggplot2::element_blank(),
                            panel.grid.minor.y = ggplot2::element_blank(),
                            panel.grid.major=ggplot2::element_blank())))
  VCF_PLOT
}
