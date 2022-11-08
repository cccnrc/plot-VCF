HG38_FASTA_FILE  <- '/home/enrico/columbia/FASTA/hg38.fa'

CHR_NAMES <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
genomestrings <- Biostrings::readDNAStringSet(HG38_FASTA_FILE)
CHRs <- genomestrings[ CHR_NAMES ]
chr_lengths <- BiocGenerics::width(CHRs)
seqinfo <- GenomeInfoDb::Seqinfo(seqnames = CHR_NAMES, seqlengths = chr_lengths, genome = 'hg38')
HG38_SEQINFO <- seqinfo

usethis::use_data(HG38_SEQINFO, overwrite = TRUE)
