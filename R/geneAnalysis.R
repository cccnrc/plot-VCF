#' Make a gene-wise analysis
#'
#' This function generates the gene-wise plot starting from the VCF input file
#'
#' @param VCF_FILE Path to the input VCF file
#' @param FASTA_FILE Path to the input FASTA file
#' @param ASSEMBLY (optional) which assembly your VCF is (hg38/hg19)
#' @param CENTILE percentile of genes to use in the plot (default is 0.9: only genes with a number of variants >90%-centile will be plotted)
#' @param SHAPE (optional) different shape for each sample
#' @param SAMPLE (optional) samples to plot (default is all samples in VCF file)
#' @param COLOR_SAMPLE (optional) list with groups and names of sample to color in the final plot
#' @param THRESHOLD (optional) a threshold line to use as Y-axis
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @param CHR_Y (optional) extend the analysis to chrY
#' @return gene-wise summary tables and plots
#' @export
geneAnalysis <- function( VCF_FILE, FASTA_FILE = FALSE, ASSEMBLY="hg38", CENTILE = 0.9, SHAPE=FALSE, SAMPLE="ALL", COLOR_SAMPLE=FALSE, THRESHOLD=FALSE, CHR_NAMES=FALSE, CHR_Y=FALSE ){
  ### check ASSEMBLY or CHR_NAMES are specified
  if ( length(CHR_NAMES) == 1 ) {
    if ( CHR_NAMES == FALSE ) {
      if ( ASSEMBLY == 'hg38' ) {
        CHR_NAMES <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
        if ( CHR_Y != FALSE ) {
          CHR_NAMES <- c( CHR_NAMES, "chrY" )
        }
      }
      else if ( ASSEMBLY == 'GRCh37' ) {
        CHR_NAMES <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
        if ( CHR_Y != FALSE ) {
          CHR_NAMES <- c( CHR_NAMES, "Y" )
        }
      } else {
        cat(' -> please specify either default assembly (hg38/GRCh37) or CHR_NAMES to plot! (see documentation)...\n')
        stop()
      }
    }
  }
  ### check COLOR_SAMPLE if passed
  if ( length(COLOR_SAMPLE) > 1 ) {
    SAMPLE_DB <- adapt_color_sample( COLOR_SAMPLE, SAMPLE )
    SAMPLE <- SAMPLE_DB$SAM
    SAMPLE_GROUP <- SAMPLE_DB$GROUP
  }
  ### output passed arguments
  cat('\n')
  cat('---------------------------------------------------------------------------\n')
  cat('------------------------>   plotVCF v0.0.0.9001   <------------------------\n')
  cat('---------------------------------------------------------------------------\n')
  cat('    -> input VCF file:', '\t', VCF_FILE, '\n' )
  cat('    -> input assembly:', '\t', ASSEMBLY, '\n' )
  cat('    -> input chromos:', '\t', CHR_NAMES, '\n' )
  cat('    -> input centile:', '\t', CENTILE, '\n' )
  cat('    -> input SAMPLE:', '\t', SAMPLE, '\n' )
  if ( FASTA_FILE != FALSE ) {
    cat('    -> input FASTA:', '\t', FASTA_FILE, '\n' )
  }
  if ( THRESHOLD != FALSE ) {
    cat('    -> input thresh:', '\t', THRESHOLD, '\n' )
  }
  if ( SHAPE != FALSE ) {
    cat('    -> input shape:', '\t', SHAPE, '\n' )
  }
  if ( CHR_Y != FALSE ) {
    cat('    -> input chrY:', '\t', CHR_Y, '\n' )
  }
  if ( length(COLOR_SAMPLE) > 0 ) {
    cat('    -> COLOR_SAMPLE:', '\t', length(COLOR_SAMPLE), 'groups specified\n' )
  }
  cat('---------------------------------------------------------------------------\n')
  cat('\n')
  ### if user specified their own FASTA file use this
  if ( FASTA_FILE != FALSE ) {
    if ( VERBOSE == TRUE ) {
      cat(' -> loading FASTA file ...\n')
    }
    SEQ <- load_chr(FASTA_FILE, CHR_NAMES=CHR_NAMES)
  } else {
    if ( ASSEMBLY == 'hg38' ) {
      for ( CHR_NAME in CHR_NAMES )
      {
        if ( ! CHR_NAME %in% names(HG38_SEQINFO) ) {
          cat(' -> ERROR: CHR_NAMES passed are not in default seqinfo for ', ASSEMBLY , ' ...\n')
          stop()
        }
      }
      SEQ <- HG38_SEQINFO[ CHR_NAMES ]
    } else if ( ASSEMBLY == 'GRCh37' ) {
      for ( CHR_NAME in CHR_NAMES )
      {
        if ( ! CHR_NAME %in% names(HG37_SEQINFO) ) {
          cat(' -> ERROR: CHR_NAMES passed are not in default seqinfo for ', ASSEMBLY , ' ...\n')
          stop()
        }
      }
      SEQ <- HG37_SEQINFO[ CHR_NAMES ]
    } else {
      cat(' -> ERROR: at least one default assembly (hg38/GRCh37) or a FASTA_FILE must be specified ...\n')
      stop()
    }
  }
  ### read the VCF file
  VCF_DATA <- read_vcf( VCF_FILE, SEQ, ASSEMBLY=ASSEMBLY, SAMPLE=SAMPLE, CHR_NAMES=CHR_NAMES )
  ### modify loaded data
  if ( (exists("SAMPLE_DB")) && ( is.data.frame(get('SAMPLE_DB')) ) ) {
    SAMPLE_DB <- SAMPLE_DB
  } else {
    SAMPLE_DB <- FALSE
  }
  RETURN_LIST <- spike_analysis( VCF_DATA, CENTILE=CENTILE, SAMPLE_DB=SAMPLE_DB, SHAPE=SHAPE, THRESHOLD=THRESHOLD, CHR_NAMES=CHR_NAMES)
  cat('\n')
  cat('---------------------------------------------------------------------------\n')
  return( RETURN_LIST )
}
