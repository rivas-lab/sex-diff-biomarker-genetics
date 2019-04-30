# filter_snps.R
# E Flynn
# 9/12/2017
#
# Code for filtering SNPs. 

FILTER_DIR = '/oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering'
OUT_DIR = '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data'

filt.file <- read.delim(sprintf("%s/variant_filter_table.tsv.gz", FILTER_DIR))
# description in file.readme - briefly filters based on HWE, MAF, manual plots
# skips displaying consequence field
#head(filt.file[,c(1:12,14:30)])

#should only include those with all_filters=0
# all filters imposes MAF of 1%, HWE < 1*10^-7, MCPI pass (if looked at)
rem.snps <- filt.file[filt.file$all_filters==0,] # 655,654 out of 784,256

# filter for ld - part of LD-pruned set
rem.snps2 <- rem.snps[rem.snps$ld_indep,] # 361,424

write.table(rem.snps2, sprintf("%s/snp_filt_metadata.txt", OUT_DIR))

write.table(rem.snps2$ID, sprintf("%s/snp_filt_list.txt", OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)