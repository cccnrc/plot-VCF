BiocManager::install()
BiocManager::install('AnnotationHub')
library(AnnotationHub)
## Load the annotation resource.
ah <- AnnotationHub()


###### GRCh38
## Query for available H.Sapiens EnsDb databases
ahDb38 <- query(ah, pattern = c("Homo Sapiens", "EnsDb"))
### latest version is 107
ens107 <- ahDb38[["AH104864"]]
seqlevelsStyle(ens107) <- "UCSC"
genes38 <- genes( ens107, columns="symbol" )
### keep only main chromosomes
CHR_NAMES <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
genes38 <- genes38[ seqnames(genes38) %in% CHR_NAMES ]
### remove missing
GENES38 <- genes38[ genes38$symbol != "" ]
genome(GENES38) <- 'hg38'
seqlevels(GENES38) <- seqlevelsInUse(GENES38)
### save object
usethis::use_data( GENES38, overwrite = TRUE )









###### GRCh37
ahDb37 <- query(ah, pattern = c("Homo Sapiens", "GRCh37"))
ens37 <- ahDb37[["AH75190"]]$codingGenome
seqlevelsStyle(ens37) <- "NCBI"
#genes37 <- genes( ens37, columns="symbol" )
### keep only main chromosomes
CHR_NAMES <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","M")
#genes38 <- genes38[ seqnames(genes38) %in% CHR_NAMES ]
### remove missing
#GENES38 <- genes38[ genes38$symbol != "" ]
### extract genes
#GENES <- c( "HLA-A", "HLA-B" )
#genes_focus <- genes38[ genes38$symbol %in% GENES ]
#names( genes_focus ) <- genes_focus$symbol
### save object
# usethis::use_data( GENES38, overwrite = TRUE )
