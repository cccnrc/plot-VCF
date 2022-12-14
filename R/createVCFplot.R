#' Make a Manhattan plot from a VCF file
#'
#' This function generates the plot starting from the VCF input file
#'
#' @param VCF_FILE Path to the input VCF file
#' @param FASTA_FILE Path to the input FASTA file
#' @param ASSEMBLY (optional) which assembly your VCF is (hg38/hg19)
#' @param GENE (optional) genes to focus the plot on
#' @param EXON (optional) plot exons
#' @param SAMPLE (optional) samples to plot (default is all samples in VCF file)
#' @param COLOR_SAMPLE (optional) list with groups and names of sample to color in the final plot
#' @param XLIM (optional) limits for Y-axis
#' @param VAR_FLAG (optional) the VCF variable to use as Y-axis for the variants, default is just position
#' @param SHAPE (optional) different shape for each sample
#' @param THRESHOLD (optional) a threshold line to use as Y-axis
#' @param ORDERED (optional) if user want to have ordered variant height in the plot or random (if value NOT specified)
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @param CHR_Y (optional) extend the analysis to chrY
#' @param SPACELINE (optional) show line between chromosomes
#' @param VERBOSE (optional) if you want steps printed to stdout
#' @return the plot
#' @export
createVCFplot <- function(VCF_FILE, FASTA_FILE = FALSE, ASSEMBLY="hg38", VAR_FLAG="POS", SHAPE=FALSE, GENE=FALSE, EXON=FALSE, SAMPLE="ALL", COLOR_SAMPLE=FALSE, XLIM=FALSE, THRESHOLD=FALSE, ORDERED=FALSE, CHR_NAMES=FALSE, CHR_Y=FALSE, SPACELINE=FALSE, VERBOSE=TRUE){
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
  ### check GENE if passed
  if ( !is.logical(get('GENE')) ) {
    removeOVERLAP <- function( GRO ) {
      GROR <- reduce(GRO)
      GRO_DF <- data.frame(GRO)
      GROR_DF <- data.frame(GROR)
      GRO_S <- GRO_DF[ GRO_DF$start %in% GROR_DF$start & GRO_DF$end %in% GROR_DF$end, ]
      GRO_S <- GRO_S[!duplicated(GRO_S$start, GRO_S$end),]
      return(makeGRangesFromDataFrame( GRO_S, keep.extra.columns=T))
    }
    ### extract genes GRanges
    GENE_GR <- GENES38[ GENES38$symbol %in% GENE ]
    names( GENE_GR ) <- GENE_GR$symbol
    seqlevels( GENE_GR ) <- seqlevelsInUse(GENE_GR)
    IDENTIFIED_GENES <- unique(names(GENE_GR))
    ### check genes are all on same chromosome
    if ( length(seqlevels( GENE_GR )) > 1 ) {
      stop("      -> all genes must be on same chromosome ...")
    }
    ### if asked for EXON analysis check only one gene was passed
    if ( is.logical(get('EXON')) ) {
      if ( EXON == TRUE ) {
        if ( length(IDENTIFIED_GENES) > 1 ) {
          stop("      -> EXON analysis can be performed only on a single gene ...")
        }
      }
    }
  }
  ### output passed arguments
  cat('\n')
  cat('---------------------------------------------------------------------------\n')
  cat('------------------------>   plotVCF v0.0.0.9001   <------------------------\n')
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
  if ( !is.logical(get('GENE')) ) {
    cat('    -> input gene(s):', '\t', GENE, '\n' )
    cat('    -> found gene(s):', '\t', IDENTIFIED_GENES, '\n' )
    if ( get('EXON') == TRUE ) {
      cat('    -> exon plot:', '\t', EXON, '\n' )
    }
  }
  if ( length(COLOR_SAMPLE) > 0 ) {
    cat('    -> COLOR_SAMPLE:', '\t', length(COLOR_SAMPLE), 'groups specified\n' )
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
  ### read the VCF file
  VCF_DATA <- read_vcf( VCF_FILE, SEQ, ASSEMBLY, SAMPLE, CHR_NAMES )
  ### modify loaded data
  LOAD_VCF_LIST <- load_vcf( VCF_DATA, SEQ, SAMPLE=SAMPLE, VAR_FLAG=VAR_FLAG )
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
  }
  if ( (exists("SAMPLE_DB")) && ( is.data.frame(get('SAMPLE_DB')) ) ) {
    SAMPLE_DB <- SAMPLE_DB
  } else {
    SAMPLE_DB <- FALSE
  }
  if ( (exists("GENE_GR")) && ( length(get('GENE_GR'))>0 ) && ( !is.logical(get('GENE_GR'))) ) {
    GENE_GR <- GENE_GR
  } else {
    GENE_GR <- FALSE
  }
  PLOT <- make_plot(VCF, SEQ, VAR_FLAG=VAR_FLAG, GENE_GR=GENE_GR, EXON=EXON, SAMPLE_DB=SAMPLE_DB, THRESHOLD=THRESHOLD, VAR_Y=VAR_Y, CHR_NAMES=CHR_NAMES, XLIM=XLIM, SHAPE=SHAPE, SPACELINE=SPACELINE)
  cat('\n')
  cat('---------------------------------------------------------------------------\n')
  PLOT
}
