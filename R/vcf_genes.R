#' return a GRanges object with gene symbols annotated
#'
#' This function takes in a GRanges object and returns it with genes annotated
#'
#' @param VCF_DATA GRanges object
#' @param COLUMN (optional) the column name of GENES38 to annotate (default: "symbol")
#' @param STRING (optional) get each gene values as comma separated string (in case of multiple genes)
#' @param SINGLE (optional) get only one genes per variant
#' @return GRanges annotated object
vcf_genes <- function( VCF_DATA, COLUMN="symbol", STRING=FALSE, SINGLE=FALSE ){
  splitColumnByOverlap <- function(query, subject, COLUMN="symbol", STRING=FALSE, SINGLE=FALSE, ...)
  {
      cat("  -> annotating VCF genes ...\n")
      olaps <- findOverlaps(query, subject, ...)
      f1 <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
      F1L <- splitAsList(mcols(query)[[COLUMN]][queryHits(olaps)], f1)
      ### return only a gene for each value
      if (( is.logical(SINGLE) ) & ( SINGLE != FALSE )) {
        F1L <- sapply(F1L,"[",1)
      }
      ### return a string element for each value
      if (( is.logical(STRING) ) & ( STRING != FALSE )) {
        F1L <- unstrsplit(F1L, sep=", ")
      }
      return( F1L )
  }
  ### check genome assembly
  VCF_ASSEMBLY <- as.character( genome(VCF_DATA)[1] )
  if (( VCF_ASSEMBLY == 'hg38' )||(VCF_ASSEMBLY == 'GRCh38')) {
    VCF_GENES <- splitColumnByOverlap(GENES38, VCF_DATA, COLUMN="symbol", STRING=STRING, SINGLE=SINGLE)
  } else {
    stop( '  -> implementation for', VCF_ASSEMBLY, 'still under developement...' )
  }
  ### add the column to VCF_DATA
  elementMetadata(VCF_DATA)[['GENE']] <- VCF_GENES
  ### actually run the function
  return( VCF_DATA )
}
