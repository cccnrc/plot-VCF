CHR_NAMES <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
CODING_LENGTH <- c( 8345447,5878924,4836344,3417120,3959062,4251743,3840011,2910726,3182390,3240638,4701713,4544590,1644627,2772368,2983594,3282459,4419593,1557059,4741161,1950206,860828,1644968,3129160,176352 )

HG38_CODING_CHR_DB <- data.frame( 'CHROM' = CHR_NAMES, 'CODING_LENGTH' = CODING_LENGTH )
HG38_CODING_CHR_DB$CHROM <- factor( HG38_CODING_CHR_DB$CHROM, levels = chr_order(HG38_CODING_CHR_DB$CHROM) )

usethis::use_data(HG38_CODING_CHR_DB, overwrite = TRUE)
