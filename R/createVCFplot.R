#' Make a Manhattan plot from a VCF file
#'
#' This function generates the plot starting from the VCF input file
#'
#' @param VCF_FILE Path to the input VCF file
#' @param FASTA_FILE Path to the input FASTA file
#' @param ASSEMBLY (optional) which assembly your VCF is (hg38/hg19)
#' @param SAMPLE (optional) samples to plot (default is all samples in VCF file)
#' @param XLIM (optional) limits for Y-axis
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants, default is just position
#' @param SHAPE (optional) different shape for each sample
#' @param THRESHOLD (optional) a threshold line to use as Y-axis
#' @param ORDERED (optional) if user want to have ordered variant height in the plot or random (if value NOT specified)
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @param CHR_Y (optional) extend the analysis to chrY
#' @param VERBOSE (optional) if you want steps printed to stdout
#' @return the plot
#' @export
createVCFplot <- function(VCF_FILE, FASTA_FILE = FALSE, ASSEMBLY="hg38", VAR_FLAG="POS", SHAPE=FALSE, SAMPLE="ALL", XLIM=FALSE, THRESHOLD=FALSE, ORDERED=FALSE, CHR_NAMES=FALSE, CHR_Y=FALSE, VERBOSE=TRUE){
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
  ### output passed arguments
  cat('\n')
  cat('---------------------------------------------------------------------------\n')
  cat('------------------------>   plotVCF v0.0.0.9000   <------------------------\n')
  cat('---------------------------------------------------------------------------\n')
  cat('    -> input VCF file:', '\t', VCF_FILE, '\n' )
  cat('    -> input assembly:', '\t', ASSEMBLY, '\n' )
  cat('    -> input chromos:', '\t', CHR_NAMES, '\n' )
  cat('    -> input VAR_FLAG:', '\t', VAR_FLAG, '\n' )
  cat('    -> input SAMPLE:', '\t', SAMPLE, '\n' )
  if ( FASTA_FILE != FALSE ) {
    cat('    -> input FASTA:', '\t', FASTA_FILE, '\n' )
  }
  if ( XLIM != FALSE ) {
    cat('    -> input XLIM:', '\t', XLIM, '\n' )
  }
  if ( THRESHOLD != FALSE ) {
    cat('    -> input thresh:', '\t', THRESHOLD, '\n' )
  }
  if ( ORDERED != FALSE ) {
    cat('    -> input order:', '\t', ORDERED, '\n' )
  }
  if ( SHAPE != FALSE ) {
    cat('    -> input shape:', '\t', SHAPE, '\n' )
  }
  if ( CHR_Y != FALSE ) {
    cat('    -> input chrY:', '\t', CHR_Y, '\n' )
  }
  cat('---------------------------------------------------------------------------\n')
  cat('\n')
  ### cannot use VAR_FLAG and ORDERED
  if (( VAR_FLAG != "POS" )&( ORDERED != FALSE )) {
    cat('\n')
    cat(' -> can not use both VAR_FLAG and ORDERED! (see documentation). Choose only one...\n')
    stop()
  }
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
  ### load VCF file
  if ( VERBOSE == TRUE ) {
    cat('\n')
    cat(' -> loading VCF file ...\n')
  }
  LOAD_VCF_LIST <- load_vcf(VCF_FILE, SEQ, ASSEMBLY=ASSEMBLY, SAMPLE=SAMPLE, VAR_FLAG=VAR_FLAG, CHR_NAMES=CHR_NAMES)
  VCF <- model_vcf(LOAD_VCF_LIST$VCF)
  SEQ <- LOAD_VCF_LIST$SEQINFO
  if ( VERBOSE == TRUE ) {
    cat(' -> arranging variants ...\n')
  }
  ### do not need those if user-specified Y axis
  VAR_Y <- FALSE
  if ( VAR_FLAG == "POS" ) {
    CHR_N <- get_chr_num(VCF, CHR_NAMES=CHR_NAMES)
    VAR_Y <- get_chr_plot_step(CHR_N, ORDERED=ORDERED)
  }
  if ( VERBOSE == TRUE ) {
    cat(' -> creating the plot ...\n')
    cat('\n')
  }
  PLOT <- make_plot(VCF, SEQ, VAR_FLAG=VAR_FLAG, THRESHOLD=THRESHOLD, VAR_Y=VAR_Y, CHR_NAMES=CHR_NAMES, XLIM=XLIM, SHAPE=SHAPE)
  cat('---------------------------------------------------------------------------\n')
  PLOT
}
