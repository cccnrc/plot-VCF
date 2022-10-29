#' Make a Manhattan plot from a VCF file
#'
#' This function generates the plot starting from the VCF input file
#'
#' @param VCF_FILE Path to the input VCF file
#' @param FASTA_FILE Path to the input FASTA file
#' @param ASSEMBLY (optional) which assembly your VCF is (hg38/hg19)
#' @param VALUE (optional) the VCF variable to use as Y-axis for the variants, default is just position
#' @param ORDERED (optional) if user want to have ordered variant height in the plot or random (if value NOT specified)
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @param VERBOSE (optional) if you want steps printed to stdout
#' @return the plot
#' @export
createVCFplot <- function(VCF_FILE, FASTA_FILE, ASSEMBLY="hg38", VALUE="POS", ORDERED=FALSE, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"), VERBOSE=TRUE){
  ### cannot use VALUE and ORDERED
  if (( VALUE != "POS" )&( ORDERED != FALSE )) {
    cat('\n')
    cat(' -> can not use both VALUE and ORDERED! (see documentation). Choose only one...\n')
  }
  if ( VERBOSE == TRUE ) {
    cat('\n')
    cat(' -> loading VCF file ...\n')
  }
  VCF <- load_vcf(VCF_FILE, ASSEMBLY)
  VCF <- model_vcf(VCF)
  if ( VERBOSE == TRUE ) {
    cat(' -> loading FASTA file ...\n')
  }
  SEQ <- load_chr(FASTA_FILE, CHR_NAMES)
  if ( VERBOSE == TRUE ) {
    cat(' -> arranging variants ...\n')
  }
  CHR_N <- get_chr_num(VCF, CHR_NAMES)
  VAR_Y <- get_chr_plot_step(CHR_N, ORDERED)
  if ( VERBOSE == TRUE ) {
    cat(' -> creating the plot ...\n')
    cat('\n')
  }
  PLOT <- make_plot(VCF, SEQ, VALUE, VAR_Y, CHR_NAMES)
  PLOT
}
