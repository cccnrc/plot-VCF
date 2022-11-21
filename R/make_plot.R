#' Create the plot based on modeled VCF
#'
#' This function defines the Y axis coordinate to plot each variant, based on the total variants in the VCF file
#'
#' @param MODELED_VCF the database from the VCF file (load_vcf() and model_vcf())
#' @param SEQINFO the object with chromosome characteristics: load_chr()
#' @param VAR_Y (optional) the Y-axis coordinates for the variants
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants
#' @param SPACELINE (optional) show line between chromosomes
#' @param EXON (optional) plot exons
#' @param GENE_GR (optional) GRanges object of genes to focus the plot on
#' @param THRESHOLD (optional) a threshold line to use as Y-axis
#' @param XLIM (optional) limits for Y-axis
#' @param SHAPE (optional) different shape for each sample
#' @param SAMPLE_DB (optional) 2-column dataframe with sample-group to color plot accordingly
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return The VCF variant plot
make_plot <- function(MODELED_VCF, SEQINFO, VAR_Y=FALSE, VAR_FLAG="POS", SPACELINE=FALSE, SHAPE=FALSE, GENE_GR=FALSE, EXON=FALSE, SAMPLE_DB=FALSE, THRESHOLD=FALSE, XLIM=FALSE, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")){
  removeOVERLAP <- function( GRO ) {
    GROR <- reduce(GRO)
    GRO_DF <- data.frame(GRO)
    GROR_DF <- data.frame(GROR)
    GRO_S <- GRO_DF[ GRO_DF$start %in% GROR_DF$start & GRO_DF$end %in% GROR_DF$end, ]
    GRO_S <- GRO_S[!duplicated(GRO_S$start, GRO_S$end),]
    return(makeGRangesFromDataFrame( GRO_S, keep.extra.columns=T))
  }
  ### check gene GRanges
  PLOT_LIMITS <- c( NA, NA )
  PLOT_TITLE <- ""
  if ( !is.logical(GENE_GR) ) {
    cat( '    -> adding GENE columns to VCF plot ... \n' )
    GENE_GR <- removeOVERLAP(GENE_GR)
    ### remove variants outside regions from MODELED_VCF
    GENE_CHROM_POS_DF <- data.frame( 'CHROM' = seqnames( GENE_GR ), 'START' = as.numeric(data.frame( ranges(GENE_GR) )$start), 'STOP' = as.numeric(data.frame( ranges(GENE_GR) )$end), 'GENE' = GENE_GR$symbol )
    GENE_CHROM_POS_DF$CHROM <- droplevels(GENE_CHROM_POS_DF$CHROM)
    MODELED_VCF$CHROM <- droplevels(MODELED_VCF$CHROM)
    MODELED_VCF_SUB <- MODELED_VCF[ as.character(MODELED_VCF$CHROM) == as.character(GENE_CHROM_POS_DF$CHROM[1]), ]
    cat( '      -> checking sample variants on specified gene(s) ... \n' )
    VECTOR_GENE <- vector()
    MODELED_VCF_SUB_EMPTY <- MODELED_VCF_SUB[0,]
    for ( i in 1:nrow(GENE_CHROM_POS_DF) )
    {
      GENE_CHROM_POS_DF_START <- GENE_CHROM_POS_DF$START[i]
      GENE_CHROM_POS_DF_STOP <- GENE_CHROM_POS_DF$STOP[i]
      GENE_CHROM_POS_DF_GENE <- GENE_CHROM_POS_DF$GENE[i]
      for ( x in 1:nrow(MODELED_VCF_SUB) )
      {
        if ((MODELED_VCF_SUB$POS[x] >= GENE_CHROM_POS_DF_START) & (MODELED_VCF_SUB$POS[x] <= GENE_CHROM_POS_DF_STOP)) {
          VECTOR_GENE <- c( VECTOR_GENE, GENE_CHROM_POS_DF_GENE )
          MODELED_VCF_SUB_EMPTY <- rbind( MODELED_VCF_SUB_EMPTY, MODELED_VCF_SUB[x,])
        }
      }
    }
    MODELED_VCF_SUB_EMPTY$GENE <- VECTOR_GENE
    MODELED_VCF_SUB <- MODELED_VCF_SUB_EMPTY
    if ( nrow(MODELED_VCF_SUB) < 1 ) {
      stop("      -> no variant found on the specified gene(s) ...")
    }
    MODELED_VCF <- MODELED_VCF_SUB
    ### keep only present genes
    GENE_PRESENT <- unique(MODELED_VCF[!is.na(MODELED_VCF$GENE),'GENE'])
    GENE_GR <- GENE_GR[ GENE_GR$symbol %in% GENE_PRESENT ]
    ### keep only valid variants in those genes
    MODELED_VCF <- MODELED_VCF[ !is.na(MODELED_VCF[,VAR_FLAG]), ]
    MODELED_VCF <- MODELED_VCF[ !is.na(MODELED_VCF[,'GENE']), ]
    ### extract gene positions
    GENE_START <- data.frame(ranges(GENE_GR))$start
    GENE_STOP <- data.frame(ranges(GENE_GR))$end
    GENE_VALUES <- c( GENE_START, GENE_STOP )
    LIMIT <- (max(GENE_VALUES)-min(GENE_VALUES))/900
    PLOT_LIMITS <- c( min(GENE_VALUES)-(LIMIT*2), max(GENE_VALUES)+(LIMIT*2) )
    ### create vectorwith start-end in succession
    GENE_LIMIT <- vector()
    for ( i in 1:length(GENE_START) )
    {
      GENE_LIMIT <- c(GENE_LIMIT, GENE_START[i])
      GENE_LIMIT <- c(GENE_LIMIT, GENE_STOP[i])
    }
    ### if asked for EXON analysis
    if ( EXON == TRUE ) {
      cat( '      -> adding EXON characteristics to VCF plot ... \n' )
      ### extract exons from saved DB
      GENE_EXONS_GR <- EXONS38[GENE_PRESENT]
      PLOT_TITLE <- paste( GENE_PRESENT[1], 'exons' )
      #GENE_EXONS_GR <- reduce(GENE_EXONS_GR)
      ### extract gene chromosome
      GENE_CHROM <- as.character(unique( MODELED_VCF[ !is.na(MODELED_VCF$CHROM), 'CHROM' ] )[1])
      ### keep only exons in main chromosome
      GENE_EXONS_GR <- GENE_EXONS_GR[ seqnames(GENE_EXONS_GR) == GENE_CHROM ]
      GENE_EXONS_GR <- removeOVERLAP( GENE_EXONS_GR )
      ### extract exons positions
      GENE_EXONS_START <- data.frame(ranges(GENE_EXONS_GR))$start
      GENE_EXONS_STOP <- data.frame(ranges(GENE_EXONS_GR))$end
      GENE_EXONS_VALUES <- c( GENE_EXONS_START, GENE_EXONS_STOP )
      EXONS_LIMIT <- (max(GENE_EXONS_VALUES)-min(GENE_EXONS_VALUES))/900
      EXONS_PLOT_LIMITS <- c( min(GENE_EXONS_VALUES)-(LIMIT*2), max(GENE_EXONS_VALUES)+(LIMIT*2) )
      ### create vectorwith start-end in succession
      GENE_EXONS_LIMIT <- vector()
      for ( i in 1:length(GENE_EXONS_START) )
      {
        GENE_EXONS_LIMIT <- c(GENE_EXONS_LIMIT, GENE_EXONS_START[i])
        GENE_EXONS_LIMIT <- c(GENE_EXONS_LIMIT, GENE_EXONS_STOP[i])
      }
      ### replace GENE PLOT values with EXON ones
      LIMIT <- EXONS_LIMIT
      PLOT_LIMITS <- EXONS_PLOT_LIMITS
      GENE_LIMIT <- GENE_EXONS_LIMIT
      GENE_START <- GENE_EXONS_START
      GENE_STOP <- GENE_EXONS_STOP
      GENE_GR <- GENE_EXONS_GR
    }
  }
  ### set Y-axis limits as specified by user
  Y_MAX <- max(MODELED_VCF[,VAR_FLAG])+((max(MODELED_VCF[,VAR_FLAG])-min(MODELED_VCF[,VAR_FLAG]))/6.6)
  if ( XLIM != FALSE ) {
    if ( length(XLIM) == 1 ) {
      XLIMITS <- c(XLIM, Y_MAX)
    } else {
      XLIMITS <- XLIM
    }
  } else {
    XLIMITS <- c(NA, Y_MAX)
  }
  ### if user specified a column to use to plot data use this
  Y_AXIS_TEXT <- ggplot2::element_blank()
  if ( VAR_FLAG != "POS" ) {
    if ( VAR_FLAG %in% colnames(MODELED_VCF) ){
      if ( 'GENE' %in% colnames(MODELED_VCF) ) {
        MODELED_VCF  <- MODELED_VCF[ MODELED_VCF[,'CHROM'] %in% CHR_NAMES, c('IND','CHROM','POS', VAR_FLAG, 'GENE')]
      } else {
        MODELED_VCF  <- MODELED_VCF[ MODELED_VCF[,'CHROM'] %in% CHR_NAMES, c('IND','CHROM','POS', VAR_FLAG)]
      }
      VAR_Y <- MODELED_VCF[,VAR_FLAG]
      Y_AXIS_TEXT <- ggplot2::element_text()
    }
    else {
      cat( "\n\t-> ERROR: undefined VCF flag selected with VAR_FLAG!\n" )
      stop()
    }
  }
  else {
    if ( 'GENE' %in% colnames(MODELED_VCF) ) {
      MODELED_VCF  <- MODELED_VCF[ MODELED_VCF[,'CHROM'] %in% CHR_NAMES, c('IND','CHROM','POS','GENE')]
    } else {
      MODELED_VCF  <- MODELED_VCF[ MODELED_VCF[,'CHROM'] %in% CHR_NAMES, c('IND','CHROM','POS')]
    }
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
  if ( !is.logical(GENE_GR) ) {
    ### remove group differences if exists
    if ( ! 'group' %in% colnames( MODELED_VCF ) ) {
      MODELED_VCF$group <- rep(1,nrow(MODELED_VCF))
    }
    ### adapt gene names
    if ( EXON == FALSE ) {
      GENE_NAMES <- as.character(data.frame(GENE_GR)$symbol)
    } else {
      GENE_NAMES <- as.character(data.frame(GENE_GR)$exon_id)
    }
    GENE_TEXT_X <- GENE_START+(( GENE_STOP - GENE_START ) / 2 )
    GENE_TEXT_Y_MIN <- max(VAR_Y)+((Y_MAX-max(VAR_Y))/2)
    GENE_TEXT_Y <- vector()
    GENE_TEXT_Y <- rep( c( GENE_TEXT_Y_MIN, Y_MAX ), length(GENE_START)/2 )
    if ( length(GENE_TEXT_Y) < length(GENE_START) ) {
      GENE_TEXT_Y <- c( GENE_TEXT_Y, GENE_TEXT_Y_MIN )
    }
    GENE_COLORS <- sample(viridis::turbo(length(GENE_START)*3))[1:length(GENE_START)]
    GENE_COLORS <- scales::hue_pal()(length(GENE_START))
    GENE_COLORS <- sample(scales::hue_pal()(length(GENE_START)*3))[1:length(GENE_START)]
    VCF_PLOT <- ggplot2::ggplot() +
                          ggplot2::geom_rect( ggplot2::aes( xmin = GENE_START-LIMIT, xmax = GENE_STOP+LIMIT, ymin = -Inf, ymax = Y_MAX ), fill = GENE_COLORS, alpha = 0.1 ) +
                          ggplot2::geom_point( data = MODELED_VCF, ggplot2::aes( x = POS, y = VAR_Y, color = group, shape = IND ), size = 1.5 ) +
                          ggplot2::geom_text( ggplot2::aes( x=GENE_TEXT_X, y=GENE_TEXT_Y, label=GENE_NAMES ), color = GENE_COLORS , angle=45, size=3 ) +
                          ggplot2::theme_bw() +
                          ggplot2::ylab( VAR_FLAG ) +
                          ggplot2::ggtitle( PLOT_TITLE ) +
                          ggplot2::scale_shape_manual( values = SHAPE_SCALE, name = '') +
                          ggplot2::scale_y_continuous( limits = XLIMITS ) +
                          ggplot2::scale_x_continuous( "POS", limits = PLOT_LIMITS, labels = scales::label_comma() ) +
                          ggplot2::guides( colour = LEGEND_COLOR, shape = LEGEND_SHAPE ) +
                          ggplot2::geom_hline( yintercept=Y_THRESHOLD, linetype="longdash", color = "red", lwd = Y_THRESHOLD_THICKNESS) +
                          ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                                          panel.background=ggplot2::element_blank(),
                                          panel.border=ggplot2::element_blank(),
                                          axis.line=ggplot2::element_blank(),
                                          plot.title = ggplot2::element_text(color="darkgrey", face="bold.italic", hjust = 0.5),
                                          axis.text.y = Y_AXIS_TEXT,
                                          axis.title.x=ggplot2::element_blank(),
                                          #axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1),
                                          axis.text.x = ggplot2::element_blank(),
                                          # axis.title.y = Y_AXIS_TEXT,
                                          legend.position = "top",
                                          legend.title=ggplot2::element_blank(),
                                          panel.grid.minor = ggplot2::element_blank(),
                                          panel.grid.major=ggplot2::element_blank())
    return( VCF_PLOT )
  }
  ### actually create the plot
  if ( "group" %in% colnames(MODELED_VCF) ) {
    VCF_PLOT <- suppressWarnings(suppressMessages(ggbio::plotGrandLinear(gr_geno, ggplot2::aes( y = VAR_Y, shape = IND, color = group ),
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
    VCF_PLOT <- suppressWarnings(suppressMessages(ggbio::plotGrandLinear(gr_geno, ggplot2::aes( x = end, y = VAR_Y, shape = IND  ),
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
                                        axis.line = ggplot2::element_blank(),
                                        axis.text.y=Y_AXIS_TEXT,
                                        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
                                        axis.title.x = ggplot2::element_blank(),
                                        axis.title.y=Y_AXIS_TEXT,
                                        legend.position = "top",
                                        panel.grid.minor.y = ggplot2::element_blank(),
                                        panel.grid.major=ggplot2::element_blank())))
  }
  cat('\n')
  return( VCF_PLOT )
}
