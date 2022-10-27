#' Load a FASTA file and get chromosomes characteristics
#'
#' This function loads the FASTA file and returns a Seqinfo object with length for each chromosome
#' or with only specified chromosomes, if CHR_NAMES is given
#'
#' @param FASTA_FILE Path to the input FASTA file
#' @param CHR_NAMES (optional) vector with chromosme names to plot
#' @return A Seqinfo object with length of each chromosome
load_chr <- function(FASTA_FILE, CHR_NAMES=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")){
  ### load FASTA
  genomestrings <- Biostrings::readDNAStringSet(FASTA_FILE)
  CHRs <- genomestrings[ CHR_NAMES ]
  chr_lengths <- BiocGenerics::width(CHRs)
  #chr_lengths <- width(genomestrings)
  seqinfo <- GenomeInfoDb::Seqinfo(seqnames = CHR_NAMES, seqlengths = chr_lengths)
  seqinfo
}
