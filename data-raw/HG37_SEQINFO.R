HG37_FASTA_FILE  <- '/home/enrico/columbia/FASTA/hs37d5.fa'

CHR_ORIGINAL_NAMES <- c( "1 dna:chromosome chromosome:GRCh37:1:1:249250621:1", "2 dna:chromosome chromosome:GRCh37:2:1:243199373:1", "3 dna:chromosome chromosome:GRCh37:3:1:198022430:1", "4 dna:chromosome chromosome:GRCh37:4:1:191154276:1", "5 dna:chromosome chromosome:GRCh37:5:1:180915260:1", "6 dna:chromosome chromosome:GRCh37:6:1:171115067:1", "7 dna:chromosome chromosome:GRCh37:7:1:159138663:1", "8 dna:chromosome chromosome:GRCh37:8:1:146364022:1", "9 dna:chromosome chromosome:GRCh37:9:1:141213431:1", "10 dna:chromosome chromosome:GRCh37:10:1:135534747:1", "11 dna:chromosome chromosome:GRCh37:11:1:135006516:1", "12 dna:chromosome chromosome:GRCh37:12:1:133851895:1", "13 dna:chromosome chromosome:GRCh37:13:1:115169878:1", "14 dna:chromosome chromosome:GRCh37:14:1:107349540:1", "15 dna:chromosome chromosome:GRCh37:15:1:102531392:1", "16 dna:chromosome chromosome:GRCh37:16:1:90354753:1" , "17 dna:chromosome chromosome:GRCh37:17:1:81195210:1" , "18 dna:chromosome chromosome:GRCh37:18:1:78077248:1" , "19 dna:chromosome chromosome:GRCh37:19:1:59128983:1" , "20 dna:chromosome chromosome:GRCh37:20:1:63025520:1" , "21 dna:chromosome chromosome:GRCh37:21:1:48129895:1" , "22 dna:chromosome chromosome:GRCh37:22:1:51304566:1" , "X dna:chromosome chromosome:GRCh37:X:1:155270560:1", "Y dna:chromosome chromosome:GRCh37:Y:2649521:59034049:1", "MT gi|251831106|ref|NC_012920.1| Homo sapiens mitochondrion, complete genome" )
CHR_NAMES <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")

genomestrings <- Biostrings::readDNAStringSet(HG37_FASTA_FILE)
CHRs <- genomestrings[ CHR_ORIGINAL_NAMES ]
names(CHRs) <- CHR_NAMES
chr_lengths <- BiocGenerics::width(CHRs)
seqinfo <- GenomeInfoDb::Seqinfo(seqnames = CHR_NAMES, seqlengths = chr_lengths, genome = 'GRCh37')
HG37_SEQINFO <- seqinfo

usethis::use_data(HG37_SEQINFO, overwrite = TRUE)
