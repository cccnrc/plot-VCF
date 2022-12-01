BiocManager::install()
BiocManager::install('AnnotationHub')
library(biomaRt)
library(AnnotationHub)
#library(Homo.sapiens)
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


###### EXONS
exons38 <- exonsBy( ens107, by='gene', columns='symbol' )
### keep only main chromosomes
exons38 <- exons38[ seqnames(exons38) %in% CHR_NAMES ]
### get gene symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- names(exons38)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
EXONS_FOUND <- vector()
EXONS_GENES <- vector()
for ( i in 1:length(exons38) )
{
  ENSID <- names(exons38)[i]
  if ( ENSID %in% G_list$ensembl_gene_id ) {
    SYMID <- G_list[ G_list$ensembl_gene_id == ENSID, 'hgnc_symbol' ][1]
    EXONS_FOUND <- c( EXONS_FOUND, i )
    EXONS_GENES <- c( EXONS_GENES, SYMID )
  }
}
exons38p <- exons38[ EXONS_FOUND ]
names(exons38p) <- EXONS_GENES
EXONS38 <- exons38p
usethis::use_data( EXONS38, overwrite = TRUE )






###### GRCh37
library(EnsDb.Hsapiens.v75)
ens37 <- genes( EnsDb.Hsapiens.v75, columns="symbol" )
# seqlevelsStyle(ens37) <- "NCBI"
### keep only main chromosomes
CHR_NAMES <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","M")
genes37 <- ens37[ seqnames(ens37) %in% CHR_NAMES ]
### remove missing

elementMetadata(genes37)[['SYMBOL']] <- sapply(genes37$SYMBOL,"[",1)
GENES37 <- genes37[ !is.na(genes37$symbol) ]
GENES37 <- GENES37[ GENES37$symbol != "" ]
### save object
usethis::use_data( GENES37, overwrite = TRUE )

###### EXONS

  query(ah, pattern = c("Homo Sapiens", "ensembl", "GRCh37"))
  ens37 <- ah[['AH10684']]
  exons37 <- ens37[ ens37$type == 'exon' ]
  exons37 <- exons37[ seqnames(exons37) %in% CHR_NAMES ]

### create a GRangesList object from single gene exons GRanges
genes37 <- unique( exons37$gene_name )
GR_VECTOR <- vector()
for ( GENE in genes37 )
{
  GENE_GR <- exons37[ exons37$gene_name == GENE ]
  GR_VECTOR <- c( GR_VECTOR, GENE_GR )
}
EXONS37_GRL <- GRangesList( GR_VECTOR )
names(EXONS37_GRL) <- genes37
### add a symbol column
mcols(EXONS37_GRL, level="within")[,'symbol'] <- mcols(EXONS37_GRL, level="within")[,'gene_name']
EXONS37 <- EXONS37_GRL
### save data
usethis::use_data( EXONS37, overwrite = TRUE )
