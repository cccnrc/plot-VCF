#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied: <input-file> <output-file> <FASTA-file>", call.=FALSE)
} else if (length(args)==3) {
  # default unordered
  args[4] = "FALSE"
}

GENO <- args[1]
OUT_PNG <-args[2]
FASTA <- args[3]
ORDERED_PLOT <- args[4]

### preset
#GENO <- '/home/enrico/columbia/IgAN-twins/jul2021/GERMLINE/FINAL-VCF-OCT2022/ITA1-ITA2-CUMC1-CUMC2-NZ.CGP.filtered.deNovo.vep.norm.VEP-ANNOTATED.SNP-filtered.PASS-only.CODING.CONTROL-unshared.AF0_0001.CONSEQUENCE.ALL.MOD.geno.tab'
#OUT_PNG <- '/home/enrico/columbia/IgAN-twins/jul2021/GERMLINE/PLOT-VCF/ALL-CHR.SNP.CTR.png'
#FASTA <- "/home/enrico/columbia/FASTA/hg38.fa"


### resumen
cat('\n')
cat('\n')
cat('\n')
cat('###################################\n')
cat('###########   PLOT VCF   ##########\n')
cat('###################################\n')
cat(paste('### --> input:', args[1], '\n', sep='\t'))
cat(paste('### --> output:', args[2], '\n', sep='\t'))
cat(paste('### --> FASTA:', args[3], '\n', sep='\t'))
cat(paste('### --> order:', args[4], '\n', sep='\t'))
cat('###################################\n')

### load libraries
cat('\n')
cat(' -> loading libraries ...\n')
cat('\n')
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(Biostrings))
suppressMessages(library(ggbio))

### remove rows where no samples has GT data
geno <- suppressWarnings(read_delim( GENO, delim = "\t", na = c("NA", "./.", '.'), comment = "#", show_col_types = FALSE ))
geno <- geno[ !is.na(geno$GT) , ]

### load FASTA
cat(' -> loading FASTA file ...\n')
cat('\n')
genomestrings <- readDNAStringSet(FASTA)
CHRs <- genomestrings[ c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY") ]
chr_lengths <- width(CHRs)
seqnames <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
#chr_lengths <- width(genomestrings)
seqinfo <- Seqinfo(seqnames = seqnames, seqlengths = chr_lengths)


### get chromosomes
cat(' -> variant number for each chromosome:\n')
CHR_N <- vector()
for ( CHR in seqnames )
{
  N <- nrow( geno[ geno$CHROM == CHR, ] )
  cat('    - ', CHR, ': ', N, '\n'  )
  CHR_N <- c( CHR_N, N )
}
cat('\n')

### get MAX value and STEP for each point
CHR_N_MAX <- max( CHR_N )
STEP <- 1/CHR_N_MAX

### assign each variant its step based on the number in the CHR
VAR_HEIGHT <- vector()
for ( TOT in CHR_N )
{
  if ( TOT > 0 ) {
    for ( i in 1:TOT )
    {
      HEIGHT_I <- STEP * i
      VAR_HEIGHT <- c( VAR_HEIGHT, HEIGHT_I )
    }
  }
}


geno2  <- geno[,c('IND','CHROM','POS')]
geno2 <- geno2 %>% arrange(factor( CHROM, levels = seqnames ))

### for unordered plot
GL <- runif(n=nrow(geno2), min=1e-12, max=.99)
# geno2 <- cbind( geno2, VAR_HEIGHT )
# head(geno2)



gr_geno <- makeGRangesFromDataFrame(geno2, keep.extra.columns = TRUE,
                                    ignore.strand = TRUE, seqinfo = seqinfo,
                                    seqnames.field = "CHROM", start.field = "POS",
                                    end.field = "POS")


### create a ordered plot or not based on user input
Y_AXIS <- VAR_HEIGHT
if ( ORDERED_PLOT != "TRUE" ) {
  if ( ORDERED_PLOT != "T" ) {
    Y_AXIS <- GL
  }
}
COL <- viridis::viridis(24)

VCF_PLOT <- suppressMessages(plotGrandLinear(gr_geno, aes(y = Y_AXIS),
                    space.skip = 0.05,
                    xlab = "Chromosome",
                    ylab = "",
                    color = COL) +
                    theme_bw() +
                    # scale_x_discrete(labels = seqinfo) +
                    theme(axis.ticks=element_blank(),
                          panel.background=element_blank(),
                          panel.border=element_blank(),
                          axis.line=element_blank(),
                          axis.text.y=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          legend.position = "none",
                          panel.grid.minor=element_blank(),
                          panel.grid.major=element_blank()))


### print the plot
png( filename = OUT_PNG,  width = 5000, height = 2500, res = 300 )
VCF_PLOT
whatever <- dev.off()






























### ENDc
