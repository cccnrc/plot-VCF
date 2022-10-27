#' Load a VCF file
#'
#' This function loads a VCF file as a matrix.
#'
#' @param VCF_FILE Path to the input VCF file
#' @return A matrix of the infile
load_vcf <- function(VCF_FILE){
  VCF_HEAD<-readLines(VCF_FILE)
  VCF_DATA<-read.table( VCF_FILE, stringsAsFactors = FALSE)

  TO_MATCH <- c("^#", "*Flag")
  MATCH <- VCF_HEAD[ grep(paste(TO_MATCH, collapse="."), VCF_HEAD ) ]
  FLAGS <- gsub('^.*ID=\\s*|\\s*,.*$', '', MATCH)

  ### replace in VCF data
  OLD_COL <- VCF_DATA[,8]
  for ( F in FLAGS ) {
    REG <- paste( F, ";", sep = '' )
    OLD_COL <- sub( REG, "", OLD_COL )
  }
  VCF_DATA[,8] <- OLD_COL

  ### column names
  VCF_HEAD <- VCF_HEAD[-(grep("#CHROM",VCF_HEAD)+1):-(length(VCF_HEAD))]
  VCF_COLNAMES <- unlist(strsplit(VCF_HEAD[length(VCF_HEAD)],"\t"))
  names(VCF_DATA) <- VCF_COLNAMES
  VCF_DATA
}
