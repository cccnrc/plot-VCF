#' Extract variant type (snp/indel)
#'
#' This function takes in a CollapsedVCF object and returns a vector: variant type (snp/indel)
#'
#' @param VCF_DATA CollapsedVCF object (as extracted by VariantAnnotation::readVcf())
#' @return non null GT vector
snp_indel <- function( VCF_DATA ){
  cat( "    -> extracting variant type ...\n")
  VCF_ISSNV <- VariantAnnotation::isSNV( VCF_DATA, singleAltOnly = FALSE )
  VCF_TYPE <- ifelse( VCF_ISSNV, "SNV", "INDEL")
  return( VCF_TYPE )
}
