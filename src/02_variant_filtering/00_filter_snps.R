# filter_snps.R
# E Flynn
# 9/12/2017
#
# Code for filtering SNPs. 

FILTER_DIR = 'variant_filtering/'
OUT_DIR = 'data/' 

require('data.table')

filt.file <- fread(sprintf("%s/variant_filter_table.tsv.gz", FILTER_DIR), data.table=FALSE)
# description in file.readme - briefly filters based on HWE, MAF, manual plots
# skips displaying consequence field
#head(filt.file[,c(1:12,14:30)])

#should only include those with all_filters=0
# all filters imposes MAF of 1%, HWE < 1*10^-7, MCPI pass (if looked at)
rem.snps <- filt.file[filt.file$all_filters==0,] # 655,654 out of 784,256

# filter for ld - part of LD-pruned set
rem.snps2 <- rem.snps[rem.snps$ld_indep,] # 361,424

# updated LD filtering later
#keep.snps <- (intersect(rem.snps$ID,prune_snps$V1))

write.table(rem.snps2, sprintf("%s/snp_filt_metadata.txt", OUT_DIR))

write.table(rem.snps2$ID, sprintf("%s/snp_filt_list.txt", OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(keep.snps, sprintf("%s/snp_filt_list2.txt", OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)