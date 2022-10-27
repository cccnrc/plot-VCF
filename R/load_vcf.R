#' Load a VCF file
#'
#' This function loads a VCF file as a matrix.
#'
#' @param VCF_FILE Path to the input VCF file
#' @param ASSEMBLY (optional) assembly of the VCF file (default hg38)
#' @return A matrix of the infile
load_vcf <- function(VCF_FILE, ASSEMBLY="hg38"){
  VCF_HEAD <- readLines(VCF_FILE)
  VCF_DATA <- suppressWarnings(VariantAnnotation::readVcf( VCF_FILE, ASSEMBLY ))

  POS_COL <- data.frame(VariantAnnotation::rowRanges(VCF_DATA)[,"paramRangeID"])[,c('seqnames', 'start')]
  VAR_COL <- data.frame(VariantAnnotation::rowRanges(VCF_DATA))[,c('QUAL', 'FILTER')]
  GENO_COL <- VariantAnnotation::geno(VCF_DATA)$GT
  VCF_BODY <- cbind( POS_COL, rownames(GENO_COL), VAR_COL, GENO_COL )
  colnames( VCF_BODY ) <- c( c( 'CHROM', 'POS', 'GT' ), colnames(VCF_BODY)[4:length(colnames(VCF_BODY))] )

  # VCF_DATA<-read.table( VCF_FILE, stringsAsFactors = FALSE)

  # TO_MATCH <- c("^#", "*Flag")
  # MATCH <- VCF_HEAD[ grep(paste(TO_MATCH, collapse="."), VCF_HEAD ) ]
  # FLAGS <- gsub('^.*ID=\\s*|\\s*,.*$', '', MATCH)

  ### replace in VCF data
  # OLD_COL <- VCF_DATA[,8]
  # for ( F in FLAGS ) {
  #   REG <- paste( F, ";", sep = '' )
  #   OLD_COL <- sub( REG, "", OLD_COL )
  # }
  # VCF_DATA[,8] <- OLD_COL

  ### column names
  # VCF_HEAD <- VCF_HEAD[-(grep("#CHROM",VCF_HEAD)+1):-(length(VCF_HEAD))]
  # VCF_COLNAMES <- unlist(strsplit(VCF_HEAD[length(VCF_HEAD)],"\t"))
  # names(VCF_DATA) <- VCF_COLNAMES
  # VCF_DATA
  VCF_BODY
}
