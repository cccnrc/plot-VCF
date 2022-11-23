#' Load a VCF file
#'
#' This function loads a VCF file as a matrix.
#'
#' @param VCF_DATA loaded VCF file (read_vcf() output)
#' @param CENTILE percentile of genes to use in the plot (default is 0.9: only genes with a number of variants >90%-centile will be plotted)
#' @param SAMPLE_DB (optional) 2-column dataframe with sample-group to color plot accordingly
#' @param SHAPE (optional) different shape for each sample - to implement
#' @param THRESHOLD (optional) a threshold line to use as Y-axis - to implement
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return gene-wise sumamries and plots
spike_analysis <- function( VCF_DATA, CENTILE = 0.9, SAMPLE_DB=FALSE, SHAPE=FALSE, THRESHOLD=FALSE, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")  ){
  ### function to extract gene values for each sample
  get_sample_gene <- function( VCF_DATA )
  {
    ### extract VCF data
    VCF_SAMPLES <- rownames(colData( VCF_DATA ))
    VCF_CHROM <- data.frame(rowRanges(VCF_DATA))$seqnames
    VCF_GT <- geno(VCF_DATA)$GT
    ### annotate VCF genes
    #VCF_GENES <- suppressWarnings(data.frame(do.call('rbind', strsplit(as.character(data.frame(info(VCF_DATA)$CSQ)$value),'|',fixed=TRUE)))[,4])
    VCF_DATA_ANNOTATED <- vcf_genes( VCF_DATA, SINGLE=TRUE )
    VCF_GENES <- elementMetadata(VCF_DATA_ANNOTATED)[['GENE']]
    ### get gene corresponding chromosome
    CHROM_GENE_SUMMARY_DB <- data.frame( 'GENE' = VCF_GENES, 'CHROM' = VCF_CHROM )
    CHROM_GENE_SUMMARY_DB <- CHROM_GENE_SUMMARY_DB[!duplicated(CHROM_GENE_SUMMARY_DB$GENE), ]
    CHROM_GENE_SUMMARY_DB <- CHROM_GENE_SUMMARY_DB[complete.cases(CHROM_GENE_SUMMARY_DB), ]
    VCF_GENES_UNIQ <- unique(CHROM_GENE_SUMMARY_DB$GENE)
    CHROM_VECTOR <- vector()
    for (GENE in VCF_GENES_UNIQ)
    {
      GENE_CHROM <- as.character(CHROM_GENE_SUMMARY_DB[ CHROM_GENE_SUMMARY_DB$GENE == GENE, 'CHROM' ])
      CHROM_VECTOR <- c( CHROM_VECTOR, GENE_CHROM )
    }
    ### extract number of non-null variants in each gene for each sample
    GENE_SUMMARY_DB <- data.frame( 'GENE' = VCF_GENES_UNIQ, 'CHROM' = CHROM_VECTOR )
    for (SAMPLE in VCF_SAMPLES)
    {
      SAMPLE_VECTOR <- vector()
      SAMPLE_GT_COL <- VCF_GT[,SAMPLE]
      SAMPLE_DB <- cbind(SAMPLE_GT_COL, data.frame( 'GENE' = VCF_GENES ))
      SAMPLE_DB <- SAMPLE_DB[ SAMPLE_DB$SAMPLE_GT_COL != ".", ]
      ### loop over single genes
      for (GENE in VCF_GENES_UNIQ)
      {
        SAMPLE_GENE_N <- nrow( SAMPLE_DB[ SAMPLE_DB$GENE == GENE, ] )
        SAMPLE_VECTOR <- c( SAMPLE_VECTOR, SAMPLE_GENE_N )
      }
      GENE_SUMMARY_DB[,SAMPLE] <- SAMPLE_VECTOR
    }
    GENE_SUMMARY_DB[,'SUM'] <- rowSums(GENE_SUMMARY_DB[,c(-1,-2)])
    return(GENE_SUMMARY_DB)
  }
  ### function to transform summary DB
  transform_summary_db <- function( SUMMARY_DB, CENTILE=0.9, SAMPLE_DB=FALSE )
  {
    ### extract samples from passed DB
    SAMPLE_VEC <- colnames(SUMMARY_DB[,c(-1,-2,-ncol(SUMMARY_DB)) ])
    ### derive each sample GROUP based on COLOR_SAMPLE
    if ( !is.logical(SAMPLE_DB) ) {
      if ( is.data.frame(SAMPLE_DB) ) {
        SAMPLE_GROUP_DB <- SAMPLE_DB
        SAMPLE_GROUP_DB$GROUP <- factor(SAMPLE_GROUP_DB$GROUP)
      }
    } else {
      SAMPLE_GROUP_DB <- data.frame( 'SAM' = SAMPLE_VEC, 'GROUP' = rep(1,length(SAMPLE_VEC)) )
    }
    ### derive DB values
    N_START <- nrow(SUMMARY_DB)
    ### keep only values above 90% centile
    SUMMARY_DB <- SUMMARY_DB[ SUMMARY_DB$SUM > as.numeric(quantile( SUMMARY_DB$SUM, CENTILE )), ]
    N_FILT <- nrow(SUMMARY_DB)
    cat( "  -> transform_summary_db() reduced DB from n.", N_START, "to n.", N_FILT, " ( centile:", CENTILE, ")\n" )
    ### add gene lenght
    GENE_VECTOR <- vector()
    for ( i in 1:nrow(SUMMARY_DB) )
    {
      GENE_NAME <- SUMMARY_DB$GENE[i]
      GENE_RANGES <- ranges(reduce(GENES38[ GENES38$symbol == GENE_NAME, ]))
      GENE_LENGTH <- as.numeric(data.frame(GENE_RANGES)$width)[1]
      GENE_VECTOR <- c( GENE_VECTOR, GENE_LENGTH )
    }
    LENGTH <- GENE_VECTOR
    SUMMARY_DB <- cbind( SUMMARY_DB, LENGTH )
    ### divide variants number for each gene length ratio
    HUMNA_GENE_MEDIAN_LENGTH <- 24000
    SUMMARY_DB$LENGTH_RATIO <- SUMMARY_DB$LENGTH / HUMNA_GENE_MEDIAN_LENGTH
    SUMMARY_DB$SUM_RATIO <- SUMMARY_DB$SUM / SUMMARY_DB$LENGTH_RATIO
    ### get divided value for each sample creating a new DF
    SAMPLE_DB <- data.frame( 'GENE' = rep(SUMMARY_DB[ , 'GENE' ], length(SAMPLE_VEC)), 'CHROM' = rep(SUMMARY_DB[ , 'CHROM' ], length(SAMPLE_VEC)) )
    SAMPLE_GENE_VECTOR <- vector()
    SAMPLE_NAME_VECTOR <- vector()
    SAMPLE_GROUP_VECTOR <- vector()
    for (SAMPLE in SAMPLE_VEC)
    {
      SAMPLE_GROUP <- SAMPLE_GROUP_DB[ SAMPLE_GROUP_DB[,'SAM'] == SAMPLE, 'GROUP' ]
      SAMPLE_GENE_VECTOR <- c( SAMPLE_GENE_VECTOR, SUMMARY_DB[,SAMPLE]/SUMMARY_DB$LENGTH_RATIO)
      SAMPLE_NAME_VECTOR <- c( SAMPLE_NAME_VECTOR, rep( SAMPLE, nrow(SUMMARY_DB) ) )
      SAMPLE_GROUP_VECTOR <- c( SAMPLE_GROUP_VECTOR, rep( SAMPLE_GROUP, nrow(SUMMARY_DB) ) )
    }
    SAMPLE_DB[ , 'RATIO' ] <- SAMPLE_GENE_VECTOR
    SAMPLE_DB[ , 'SAMPLE' ] <- SAMPLE_NAME_VECTOR
    SAMPLE_DB[ , 'GROUP' ] <- SAMPLE_GROUP_VECTOR
    ### order returned database
    SUMMARY_DB[ , 'GENE' ] <- factor(SUMMARY_DB[ , 'GENE' ], levels = SUMMARY_DB[ , 'GENE' ])
    SUMMARY_DB[ , 'CHROM' ] <- factor(SUMMARY_DB[ , 'CHROM' ], levels = chr_order(SUMMARY_DB[ , 'CHROM' ]))
    SAMPLE_DB[ , 'GENE' ] <- factor(SAMPLE_DB[ , 'GENE' ], levels = levels(SUMMARY_DB[ , 'GENE' ]))
    SAMPLE_DB[ , 'CHROM' ] <- factor(SAMPLE_DB[ , 'CHROM' ], levels = levels(SUMMARY_DB[ , 'CHROM' ]))
    CHROM_LEV <- levels( SUMMARY_DB[ , 'CHROM' ] )
    CHR_VECTOR <- vector()
    for ( CHR in CHR_NAMES )
    {
      if ( CHR %in% CHROM_LEV ) {
        CHR_VECTOR <- c( CHR_VECTOR, CHR )
      }
    }
    SUMMARY_DB[ , 'CHROM' ] <- factor(SUMMARY_DB[ , 'CHROM' ], levels = chr_order(CHR_VECTOR))
    SAMPLE_DB[ , 'CHROM' ] <- factor(SAMPLE_DB[ , 'CHROM' ], levels = chr_order(CHR_VECTOR))
    ### create GROUP DB to plot out, if asked
    if ( !is.logical(SAMPLE_DB) ) {
      ### prepare the DB
      SUMMARY_GROUP_GENE_VEC <- vector()
      SUMMARY_GROUP_CHROM_VEC <- vector()
      for ( i in 1:nrow(SUMMARY_DB) )
      {
        SUMMARY_GROUP_GENE_VEC <- c( SUMMARY_GROUP_GENE_VEC, rep( as.character(SUMMARY_DB[i,'GENE']), length(levels(SAMPLE_GROUP_DB$GROUP))) )
        SUMMARY_GROUP_CHROM_VEC <- c( SUMMARY_GROUP_CHROM_VEC, rep( as.character(SUMMARY_DB[i,'CHROM']), length(levels(SAMPLE_GROUP_DB$GROUP))) )

      }
      SUMMARY_GROUP_DB <- data.frame( 'GENE' = factor(SUMMARY_GROUP_GENE_VEC), 'CHROM' = factor(SUMMARY_GROUP_CHROM_VEC) )
      SUMMARY_GROUP_DB$GENE <- factor( SUMMARY_GROUP_DB$GENE, levels = levels(SUMMARY_DB$GENE) )
      SUMMARY_GROUP_DB$CHROM <- factor( SUMMARY_GROUP_DB$CHROM, levels = levels(SUMMARY_DB$CHROM) )
      GROUP_GENE_VEC <- vector()
      GROUP_VEC <- vector()
      for ( GENE in levels(SUMMARY_GROUP_DB$GENE) )
      {
        for ( GROUP in levels(SAMPLE_GROUP_DB$GROUP) )
        {
          GROUP_VEC <- c( GROUP_VEC, GROUP )
          GROUP_N <- 0
          for ( SAMPLE in SAMPLE_GROUP_DB[ SAMPLE_GROUP_DB[,'GROUP'] == GROUP, 'SAM'] )
          {
            GROUP_N <- GROUP_N + sum( SUMMARY_DB[ SUMMARY_DB[,'GENE'] == GENE, SAMPLE ], na.rm = TRUE )
          }
          GROUP_GENE_VEC <- c( GROUP_GENE_VEC, GROUP_N / SUMMARY_DB[ SUMMARY_DB$GENE == GENE, ]$LENGTH_RATIO )
        }
      }
      SUMMARY_GROUP_DB[,'GROUP'] <- GROUP_VEC
      SUMMARY_GROUP_DB[,'RATIO'] <- GROUP_GENE_VEC
      SUMMARY_GROUP_DB[,'GROUP'] <- factor( SUMMARY_GROUP_DB[,'GROUP'] )
      ### return 3 DB
      SUMMARY_RETURN <- list( 'SUM' = SUMMARY_DB, 'RATIO' = SAMPLE_DB, 'GROUP' = SUMMARY_GROUP_DB )
      return(SUMMARY_RETURN)
    }
    ### return both DB
    SUMMARY_RETURN <- list( 'SUM' = SUMMARY_DB, 'RATIO' = SAMPLE_DB )
    return(SUMMARY_RETURN)
  }
  ### function to plot spike summaries
  make_summary_plot <- function( SUMMARY_DB_TR, SAMPLE_DB=FALSE, THRESHOLD = FALSE, SHAPE = FALSE )
  {
    ### maximum Y value for plots
    SUM_Y_MAX <- max( SUMMARY_DB_TR$SUM$SUM_RATIO )
    SAMPLE_Y_MAX <- max( SUMMARY_DB_TR$RATIO$RATIO )
    if ( !is.logical(SAMPLE_DB) ) {
      COLOR_Y_MAX <- max( SUMMARY_DB_TR$GROUP$RATIO )
    }
    ### calculate chromosome borders
    CHROM_LEVELS <- levels(SUMMARY_DB_TR$SUM$CHROM)
    CHROM_LINES <- vector()
    CHROM_TEXT <- vector()
    for ( i in 1:length(CHROM_LEVELS) )
    {
      CHROM_N <- nrow( SUMMARY_DB_TR$SUM[SUMMARY_DB_TR$SUM[,'CHROM'] == CHROM_LEVELS[i],] )
      if ( i == 1 ) {
        CHROM_LINES <- c( CHROM_LINES, CHROM_N + 0.5 )
        CHROM_TEXT <- c( CHROM_TEXT, (CHROM_N+0.5)/2 )
      } else {
        CHROM_LINES <- c( CHROM_LINES, CHROM_N + CHROM_LINES[i-1] )
        CHROM_TEXT <- c( CHROM_TEXT, (CHROM_N/2) + CHROM_LINES[i-1] )
      }
    }
    CHROM_TEXT_DB <- data.frame( 'CHROM' = factor(CHROM_LEVELS, levels = levels(SUMMARY_DB_TR$SUM$CHROM)), 'X' = as.numeric(CHROM_TEXT) )
    ### if COLOR_SAMPLE create the relative plot
    if ( !is.logical(SAMPLE_DB) ) {
      SAMPLE_GROUP_PLOT <- ggplot2::ggplot(SUMMARY_DB_TR$GROUP) +
                                  ggplot2::geom_line(ggplot2::aes(x=GENE, y=RATIO, color=GROUP, group=GROUP)) +
                                  ggplot2::geom_vline( xintercept = CHROM_LINES, linetype = "dotdash", color = "#777777", size = 0.2 ) +
                                  ggplot2::geom_text( data = CHROM_TEXT_DB, ggplot2::aes( x = X, y = 1.05*COLOR_Y_MAX, label = levels(CHROM) ),  color = "#777777", size = 3 ) +
                                  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
                                  ggplot2::xlab("") +
                                  ggplot2::ylab("variant (n) / gene length (ratio)") +
                                  ggplot2::theme_bw() +
                                  ggplot2::theme(
                                    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                                    panel.grid.major.x = ggplot2::element_blank() ,
                                    panel.grid.minor.x = ggplot2::element_blank() ,
                                    panel.grid.minor.y = ggplot2::element_blank(),
                                    legend.position="top",
                                    legend.title = ggplot2::element_blank()
                                  )
    }
    SAMPLE_PLOT <- ggplot2::ggplot(SUMMARY_DB_TR$RATIO) +
                                ggplot2::geom_line(ggplot2::aes(x=GENE, y=RATIO, color=SAMPLE, group=SAMPLE)) +
                                ggplot2::geom_vline( xintercept = CHROM_LINES, linetype = "dotdash", color = "#777777", size = 0.2 ) +
                                ggplot2::geom_text( data = CHROM_TEXT_DB, ggplot2::aes( x = X, y = 1.05*SAMPLE_Y_MAX, label = levels(CHROM) ),  color = "#777777", size = 3 ) +
                                ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
                                ggplot2::xlab("") +
                                ggplot2::ylab("variant (n) / gene length (ratio)") +
                                ggplot2::theme_bw() +
                                ggplot2::theme(
                                  axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                                  panel.grid.major.x = ggplot2::element_blank() ,
                                  panel.grid.minor.x = ggplot2::element_blank() ,
                                  panel.grid.minor.y = ggplot2::element_blank(),
                                  legend.position="top",
                                  legend.text = ggplot2::element_text(size = 4),
                                  legend.title = ggplot2::element_blank()
                                )
    CHROM_COLORS <- sample(viridis::turbo(length(levels(SUMMARY_DB_TR$SUM$CHROM))*3))[1:length(levels(SUMMARY_DB_TR$SUM$CHROM))]
    SUM_PLOT <- ggplot2::ggplot(SUMMARY_DB_TR$SUM) +
                                ggplot2::geom_line(ggplot2::aes(x=GENE, y=SUM_RATIO, color=CHROM, group=1)) +
                                ggplot2::geom_vline( xintercept = CHROM_LINES, linetype = "dotdash", color = "#777777", size = 0.2 ) +
                                ggplot2::geom_text( data = CHROM_TEXT_DB, ggplot2::aes( x = X, y = 1.05*SUM_Y_MAX, label = levels(CHROM), color = levels(CHROM) ), size = 3 ) +
                                #ggplot2::scale_color_manual( values = CHROM_COLORS ) +
                                ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
                                ggplot2::xlab("") +
                                ggplot2::ylab("variant (n) / gene length (ratio)") +
                                ggplot2::theme_bw() +
                                ggplot2::theme(
                                  axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                                  panel.grid.major.x = ggplot2::element_blank() ,
                                  panel.grid.minor.x = ggplot2::element_blank() ,
                                  panel.grid.minor.y = ggplot2::element_blank(),
                                  legend.position="none",
                                  legend.title = ggplot2::element_blank()
                                )
    if ( !is.logical(SAMPLE_DB) ) {
      MERGED_PLOT <- ggpubr::ggarrange(SUM_PLOT, SAMPLE_GROUP_PLOT, SAMPLE_PLOT,
                                        labels = c("sum", "group", "sample"),
                                        ncol = 1, nrow = 3)
      PLOT_RETURN <- list( 'SAMPLE' = SAMPLE_PLOT, 'GROUP' = SAMPLE_GROUP_PLOT, 'SUM' = SUM_PLOT, 'MERGED' = MERGED_PLOT )
    } else {
      MERGED_PLOT <- ggpubr::ggarrange(SUM_PLOT, SAMPLE_PLOT,
                                        labels = c("sum", "sample"),
                                        ncol = 1, nrow = 2)
      PLOT_RETURN <- list( 'SAMPLE' = SAMPLE_PLOT, 'SUM' = SUM_PLOT, 'MERGED' = MERGED_PLOT )
    }
  }
  ### run functions
  SUMMARY <- get_sample_gene( VCF_DATA )
  SUMMARY_TR <- transform_summary_db( SUMMARY, CENTILE=CENTILE, SAMPLE_DB=SAMPLE_DB )
  SUMMARY_PLOTS <- make_summary_plot( SUMMARY_TR, SAMPLE_DB=SAMPLE_DB, THRESHOLD=THRESHOLD, SHAPE=SHAPE )
  names(SUMMARY_TR)[2] <- "SAMPLE"
  RETURN_LIST <- list( 'TAB' = SUMMARY_TR, 'PLOT' = SUMMARY_PLOTS )
  return( RETURN_LIST )
}
