#' Chromosome analysis of passed VCF file
#'
#' This function takes in the loaded VCF (read_vcf()) and operates calculations for the chromosome analysis
#'
#' @param VCF_DATA loaded VCF file (read_vcf() output)
#' @param METHOD how to count (and plot) variants for each chromosome ("RAW"=raw number, "LEN"=n/chr-length(ratio), "COD"=n/chr-coding-length(ratio), "LENCOD"=n/chr-length(ratio)*chr-coding-length(ratio))
#' @param SAMPLE_DB (optional) 2-column dataframe with sample-group to color plot accordingly
#' @param SHAPE (optional) different shape for each sample - to implement
#' @param THRESHOLD (optional) a threshold line to use as Y-axis - to implement
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return gene-wise sumamries and plots
chr_analysis <- function( VCF_DATA, METHOD="RAW", SAMPLE_DB=FALSE, SHAPE=FALSE, THRESHOLD=FALSE, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")  ){

  ### function to extract chromosome values for each sample
  get_sample_chr <- function( VCF_DATA, METHOD="RAW" )
  {
    ### extract VCF data
    VCF_SAMPLES <- rownames(colData( VCF_DATA ))
    VCF_CHROM <- data.frame(rowRanges(VCF_DATA))$seqnames
    VCF_GT <- geno(VCF_DATA)$GT
    ### prepare chromosome DB
    CHR_DF <- data.frame( 'CHROM' = chr_order(VCF_CHROM) )
    CHR_DF$CHROM <- factor( CHR_DF$CHROM, levels = chr_order(CHR_DF$CHROM) )
    ### loop over sample to get variants in each chromosome
    for ( SAMPLE in VCF_SAMPLES )
    {
      ### recreate a 2-cols DB with CHR and GT
      SAMPLE_VECTOR <- vector()
      SAMPLE_DB <- data.frame( 'CHROM' = VCF_CHROM, 'GT' = VCF_GT[,SAMPLE] )
      for ( CHR in levels(CHR_DF$CHROM) )
      {
        ### extract non-null variants in that chromosome for each sample
        SAMPLE_CHR_GT <- SAMPLE_DB[ SAMPLE_DB$CHROM == CHR, 'GT' ]
        SAMPLE_CHR_GT_NONULL <- gt_null( SAMPLE_CHR_GT )
        SAMPLE_CHR_GT_NONULL_N <- length(SAMPLE_CHR_GT_NONULL)
        SAMPLE_VECTOR <- c( SAMPLE_VECTOR, SAMPLE_CHR_GT_NONULL_N )
      }
      ### add sample column to chromosome DB
      CHR_DF[,SAMPLE] <- SAMPLE_VECTOR
    }
    ### add SUM column (all samples sum)
    CHR_DF[,'SUM'] <- rowSums( CHR_DF[,-1] )

    ### based on specified METHOD define column output
    cat( '    -> passed chromosomes plotting method:', METHOD, '\n' )
    if ( METHOD == "RAW" ) {
      return( CHR_DF )
    }

    ### check genome assembly
    VCF_ASSEMBLY <- as.character( genome(VCF_DATA)[1] )
    if (( VCF_ASSEMBLY == 'hg38' )||(VCF_ASSEMBLY == 'GRCh38')) {
      ASSEMBLY_SEQINFO <- HG38_SEQINFO
    } else  if (( VCF_ASSEMBLY == 'hg19' )||(VCF_ASSEMBLY == 'GRCh37')) {
      ASSEMBLY_SEQINFO <- HG37_SEQINFO
    } else {
      stop( '  -> implementation for', VCF_ASSEMBLY, 'still under developement...' )
    }
    ### extract each chromosome length and cosing-length and add to CHR_DF
    CHR_LENGTH_VECTOR <- vector()
    CHR_CODING_LENGTH_VECTOR <- vector()
    for ( CHR in levels(CHR_DF$CHROM) )
    {
      CHR_LENGTH <- seqlengths(ASSEMBLY_SEQINFO)[CHR]
      CHR_CODING_LENGTH <- HG38_CODING_CHR_DB[ HG38_CODING_CHR_DB$CHROM == CHR, 'CODING_LENGTH' ]
      CHR_LENGTH_VECTOR <- c( CHR_LENGTH_VECTOR, CHR_LENGTH )
      CHR_CODING_LENGTH_VECTOR <- c( CHR_CODING_LENGTH_VECTOR, CHR_CODING_LENGTH )
    }
    CHR_DF[,'LENGTH'] <- CHR_LENGTH_VECTOR
    CHR_DF[,'CODING_LENGTH'] <- CHR_CODING_LENGTH_VECTOR
    CHR_DF[,'CODING_RATIO'] <- CHR_DF[,'CODING_LENGTH'] / CHR_DF[,'LENGTH']
    CHR_DF[,'CODING_RATE'] <- CHR_DF[,'CODING_LENGTH'] / median(CHR_DF[,'CODING_LENGTH'])
    ### extract median chromosome length
    MEDIAN_CHR_LENGTH <- median( seqlengths(ASSEMBLY_SEQINFO)[levels(CHR_DF$CHROM)] )
    cat( '    -> median length of passed chromosomes:', MEDIAN_CHR_LENGTH, 'bp\n' )
    ### calculate each chromosome length and coding-length ratio
    CHR_DF[,'LENGTH_RATIO'] <- CHR_DF$LENGTH / MEDIAN_CHR_LENGTH
    CHR_DF[,'CODING_LENGTH_RATIO'] <- CHR_DF[,'CODING_RATE'] / CHR_DF[,'LENGTH_RATIO']

    ### define output based on METHOD
    if ( METHOD == "LEN" ) {
      for ( SAMPLE in VCF_SAMPLES )
      {
        CHR_DF[,SAMPLE] <- CHR_DF[,SAMPLE] / CHR_DF[,'LENGTH_RATIO']
      }
      CHR_DF[,'SUM'] <- CHR_DF[,'SUM'] / CHR_DF[,'LENGTH_RATIO']
      CHR_DF <- CHR_DF[ ,  -which(names(CHR_DF) %in% c( 'LENGTH', 'CODING_LENGTH', 'CODING_RATIO', 'CODING_RATE', 'LENGTH_RATIO', 'CODING_LENGTH_RATIO' )) ]
    } else if ( METHOD == "COD" ) {
      for ( SAMPLE in VCF_SAMPLES )
      {
        CHR_DF[,SAMPLE] <- CHR_DF[,SAMPLE] / CHR_DF[,'CODING_RATE']
      }
      CHR_DF[,'SUM'] <- CHR_DF[,'SUM'] / CHR_DF[,'CODING_RATE']
      CHR_DF <- CHR_DF[ , -which(names(CHR_DF) %in% c( 'LENGTH', 'CODING_LENGTH', 'CODING_RATIO', 'CODING_RATE', 'LENGTH_RATIO', 'CODING_LENGTH_RATIO' )) ]
    } else if ( METHOD == "LENCOD" ) {
      for ( SAMPLE in VCF_SAMPLES )
      {
        CHR_DF[,SAMPLE] <- CHR_DF[,SAMPLE] * CHR_DF[,'CODING_LENGTH_RATIO']
      }
      CHR_DF[,'SUM'] <- CHR_DF[,'SUM'] * CHR_DF[,'CODING_LENGTH_RATIO']
      CHR_DF <- CHR_DF[ ,  -which(names(CHR_DF) %in% c( 'LENGTH', 'CODING_LENGTH', 'CODING_RATIO', 'CODING_RATE', 'LENGTH_RATIO', 'CODING_LENGTH_RATIO' )) ]
    }
    ### annotate VCF genes
    return( CHR_DF )
  }


  ### create DB as needed by ggplot2
  transform_summary_chr_db <- function(SUMMARY_DB, SAMPLE_DB=FALSE) {
    ### extract overall chromosome DB
    SUMMARY_DB_SUM <- SUMMARY_DB[,c('CHROM','SUM')]
    ### extract single-sample chromosome DB
    SAMPLE_VECTOR <- vector()
    CHROM_VECTOR <- vector()
    VALUE_VECTOR <- vector()
    for ( i in 2:(ncol(SUMMARY_DB)-1) )
    {
      SAMPLE <- as.character(colnames(SUMMARY_DB)[i])
      SAMPLE_VECTOR <- c( SAMPLE_VECTOR, rep( SAMPLE, nrow(SUMMARY_DB) ) )
      CHROM_VECTOR <- c( CHROM_VECTOR, as.character(SUMMARY_DB[,'CHROM']) )
      for ( r in 1:nrow(SUMMARY_DB) )
      {
        SAMPLE_CHR_N <- SUMMARY_DB[ r, SAMPLE ]
        VALUE_VECTOR <- c( VALUE_VECTOR, SAMPLE_CHR_N )
      }
    }
    SUMMARY_DB_SINGLE <- data.frame( 'CHROM' = CHROM_VECTOR, 'SAMPLE' = SAMPLE_VECTOR, 'VALUE' = VALUE_VECTOR )
    SUMMARY_DB_SINGLE$CHROM <- factor(SUMMARY_DB_SINGLE$CHROM, levels = chr_order(SUMMARY_DB_SINGLE$CHROM))

    ### create GROUP DB to plot out, if asked
    if ( !is.logical(SAMPLE_DB) ) {
      SAMPLE_DB$GROUP <- factor(SAMPLE_DB$GROUP)
      ### prepare the DB: CHROM column
      SUMMARY_GROUP_CHROM_VEC <- vector()
      for ( i in 1:nrow(SUMMARY_DB) )
      {
        SUMMARY_GROUP_CHROM_VEC <- c( SUMMARY_GROUP_CHROM_VEC, rep( as.character(SUMMARY_DB[i,'CHROM']), length(levels(SAMPLE_DB$GROUP))) )
      }
      SUMMARY_GROUP_DB <- data.frame( 'CHROM' = factor(SUMMARY_GROUP_CHROM_VEC) )
      SUMMARY_GROUP_DB$CHROM <- factor( SUMMARY_GROUP_DB$CHROM, levels = chr_order(SUMMARY_GROUP_DB$CHROM) )
      ### extract values for each group
      GROUP_CHROM_VEC <- vector()
      GROUP_VEC <- vector()
      for ( CHR in levels(SUMMARY_GROUP_DB$CHROM) )
      {
        for ( GROUP in levels(SAMPLE_DB$GROUP) )
        {
          GROUP_VEC <- c( GROUP_VEC, GROUP )
          GROUP_N <- 0
          for ( SAMPLE in SAMPLE_DB[ SAMPLE_DB[,'GROUP'] == GROUP, 'SAM'] )
          {
            GROUP_N <- GROUP_N + SUMMARY_DB[ SUMMARY_DB[,'CHROM'] == CHR, SAMPLE ]
          }
          GROUP_CHROM_VEC <- c( GROUP_CHROM_VEC, GROUP_N )
        }
      }
      SUMMARY_GROUP_DB[,'GROUP'] <- GROUP_VEC
      SUMMARY_GROUP_DB[,'VALUE'] <- GROUP_CHROM_VEC
      SUMMARY_GROUP_DB[,'GROUP'] <- factor( SUMMARY_GROUP_DB[,'GROUP'] )
      ### return 3 DB
      SUMMARY_RETURN <- list( 'SUM' = SUMMARY_DB_SUM, 'SINGLE' = SUMMARY_DB_SINGLE, 'GROUP' = SUMMARY_GROUP_DB )
      return(SUMMARY_RETURN)
    }

    ### return both DB
    SUMMARY_RETURN <- list( 'SUM' = SUMMARY_DB_SUM, 'SINGLE' = SUMMARY_DB_SINGLE )
    return(SUMMARY_RETURN)
  }



  ### function to plot chromosome summaries
  make_summary_chr_plot <- function( SUMMARY_RETURN, METHOD="RAW", SAMPLE_DB=FALSE, THRESHOLD = FALSE, SHAPE = FALSE )
  {
    ### create Y-label based on METHOD passed
    if ( METHOD == "RAW" ) {
      YLABEL <- "variant (n)"
    } else if ( METHOD == "LEN" ) {
      YLABEL <- "variant (n) / chr-length (ratio)"
    } else if ( METHOD == "COD" ) {
      YLABEL <- "variant (n) / chr-coding (ratio)"
    } else if ( METHOD == "LENCOD" ) {
      YLABEL <- "variant (n) / chr-coding-length (ratio)"
    }
    ### if COLOR_SAMPLE create the relative plot
    if ( !is.logical(SAMPLE_DB) ) {
      SAMPLE_GROUP_PLOT <- ggplot2::ggplot(SUMMARY_RETURN$GROUP) +
                                  ggplot2::geom_line(ggplot2::aes(x=CHROM, y=VALUE, color=GROUP, group=GROUP)) +
                                  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
                                  ggplot2::xlab("") +
                                  ggplot2::ylab(YLABEL) +
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
    SAMPLE_PLOT <- ggplot2::ggplot(SUMMARY_RETURN$SINGLE) +
                                ggplot2::geom_line(ggplot2::aes(x=CHROM, y=VALUE, color=SAMPLE, group=SAMPLE)) +
                                ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
                                ggplot2::xlab("") +
                                ggplot2::ylab(YLABEL) +
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
    CHROM_COLORS <- sample(viridis::turbo(length(levels(SUMMARY_RETURN$SUM$CHROM))*3))[1:length(levels(SUMMARY_RETURN$SUM$CHROM))]
    SUM_PLOT <- ggplot2::ggplot(SUMMARY_RETURN$SUM) +
                                ggplot2::geom_line(ggplot2::aes(x=CHROM, y=SUM, color=CHROM, group=1)) +
                                #ggplot2::scale_color_manual( values = CHROM_COLORS ) +
                                ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
                                ggplot2::xlab("") +
                                ggplot2::ylab(YLABEL) +
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
  SUMMARY <- get_sample_chr( VCF_DATA, METHOD=METHOD )
  SUMMARY_TR <- transform_summary_chr_db( SUMMARY, SAMPLE_DB=SAMPLE_DB )
  SUMMARY_PLOTS <- make_summary_chr_plot( SUMMARY_TR, METHOD=METHOD, SAMPLE_DB=SAMPLE_DB, THRESHOLD=THRESHOLD, SHAPE=SHAPE )
  names(SUMMARY_TR)[2] <- "SAMPLE"
  RETURN_LIST <- list( 'TAB' = SUMMARY_TR, 'PLOT' = SUMMARY_PLOTS )
  return( RETURN_LIST )
}
