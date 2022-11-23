#' Read a VCF file
#'
#' This function read a VCF file and return the readVcf() object
#'
#' @param VCF_FILE Path to the input VCF file
#' @param SEQINFO the object with chromosome characteristics: load_chr()
#' @param ASSEMBLY (optional) assembly of the VCF file (default hg38)
#' @param SAMPLE (optional) samples to plot (default is all samples in VCF file)
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return A matrix of the infile
read_vcf <- function( VCF_FILE, SEQINFO, ASSEMBLY="hg38", SAMPLE="ALL", CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY") ){
  ### prepare VCF_FILE for reading
  getExtension <- function(FILE_PATH) {
      ex <- strsplit(basename(FILE_PATH), split="\\.")[[1]]
      cat('  -> passed VCF file extension:', ex[length(ex)], '\n')
      return(ex[length(ex)])
  }
  compressVCF <- function( VCF_PATH ) {
    cat('  -> compressing VCF file ...\n')
    VCF_GZ <- Rsamtools::bgzip( VCF_PATH, tempfile() )
    return( VCF_GZ )
  }
  indexVCFGZ <- function( VCF_GZ ) {
    cat('  -> indexing VCF file ...\n')
    idx <- Rsamtools::indexTabix( VCF_GZ, "vcf" )
    tab <- Rsamtools::TabixFile( VCF_GZ, idx )
    return( tab )
  }
  getVCFchr <- function ( VCF_GZ ) {
    CHR_VCF_GZ <- Rsamtools::headerTabix( VCF_GZ, "vcf" )$seqnames
    return(CHR_VCF_GZ)
  }
  getPARAM <- function( CHR_NAMES, SEQINFO, SAMPLE ) {
    SEQ_CHR <- SEQINFO[ CHR_NAMES ]
    GR_CHR <- GenomicRanges::GRanges(SEQ_CHR)
    PARAM <- VariantAnnotation::ScanVcfParam( which = GR_CHR )
    if (( length(SAMPLE) > 1 ) || ( SAMPLE != "ALL" )) {
      PARAM <- VariantAnnotation::ScanVcfParam( which = GR_CHR, samples = SAMPLE )
    }
    return( PARAM )
  }

  ### compress and index file (if not done yet)
  VCF_EXT <- getExtension(VCF_FILE)
  if (( VCF_EXT == "vcf" ) || ( VCF_EXT == "VCF" )) {
    #cat(' -> passed file is .vcf Compressing and indexing it ...\n')
    VCF_GZ <- compressVCF( VCF_FILE )
    VCF_TAB <- indexVCFGZ( VCF_GZ )
    ### check and get chromosmes actually in VCF file
    VCF_CHR_NAMES <- getVCFchr( VCF_GZ )
  } else if (( VCF_EXT == "gz" ) || ( VCF_EXT == "GZ" )) {
    #cat(' -> passed file is .vcf.gz Indexing it ...\n')
    VCF_TAB <- indexVCFGZ( VCF_FILE )
    ### check and get chromosmes actually in VCF file
    VCF_CHR_NAMES <- getVCFchr( VCF_FILE )
  } else {
    cat('\n')
    cat(' -> passed VCF file is neither .vcf or .vcf.gz! (see documentation)...\n')
    stop()
  }
  ### compare CHR_NAMES with chromosomes actually in VCF file
  CHR_ABS <- setdiff( CHR_NAMES, VCF_CHR_NAMES )
  if ( length( CHR_ABS ) > 0 ) {
    for ( CHR_AB in CHR_ABS ) {
      cat( "    -> VCF, no variant found in:", CHR_AB,"(will remove) \n" )
    }
    CHR_NAMES <- CHR_NAMES[!CHR_NAMES %in%CHR_ABS]
  }
  ### adapt SEQINFO
  SEQINFO <- SEQINFO[ CHR_NAMES ]
  ### actually read and load VCF file
  VCF_PARAM <- getPARAM( CHR_NAMES, SEQINFO, SAMPLE )
  cat('  -> reading VCF file ...\n')
  VCF_DATA <- suppressWarnings(VariantAnnotation::readVcf( VCF_TAB, ASSEMBLY, VCF_PARAM ))
  ### return loaded VCF file
  return( VCF_DATA )
}
