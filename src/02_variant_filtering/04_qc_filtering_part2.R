# 04_qc_filtering_part2.R
# E Flynn
#
# Code for QC filtering X, XY, Y, and MT chromosomes.

require('tidyverse')
require('data.table')
DATA.DIR <- "data/"
options(stringsAsFactors=FALSE)

# read in the table
vars.to.keep <- read.table(sprintf("%s/snp_filt_list.txt" , DATA.DIR), header=FALSE)

# read in all the chr qc
cX_ld <- read.table(sprintf("%s/chr_qc/chrX_ld.prune.in", DATA.DIR), header=FALSE)
cXY_ld <- read.table(sprintf("%s/chr_qc/chrXY_ld.prune.in", DATA.DIR), header=FALSE)
cX_tab <- read.table(sprintf("%s/chr_qc/chrX_qc_table.txt", DATA.DIR), sep=" ", header=TRUE)
cXY_tab <- read.table(sprintf("%s/chr_qc/chrXY_qc_table.txt", DATA.DIR), sep=" ", header=TRUE)
cXY_hardy <- read_tsv(sprintf("%s/chr_qc/chrXY_hardy.hardy", DATA.DIR))
cX_hardy <- read_tsv(sprintf("%s/chr_qc/chrX_hardy.hardy.x", DATA.DIR))

# combine data
combinedX <- cX_tab %>% mutate(LD2=ifelse(SNP %in% cX_ld$V1,1,0))
combinedXY <- cXY_tab %>% mutate(LD2=ifelse(SNP %in% cXY_ld$V1,1,0))
combinedX2 <- combinedX %>% select(-LD, -keep) %>% rename(LD=LD2) 
combinedXY2 <- combinedXY %>% select(-LD, -keep) %>% rename(LD=LD2) 

# filter
cXY_h <- left_join(combinedXY2, cXY_hardy %>% select(ID, P), by=c("SNP"="ID")) %>% rename(HWE_P=P)
cX_h <- left_join(combinedX2, cX_hardy %>% select(ID, P), by=c("SNP"="ID")) %>% rename(HWE_P=P)
cX_keep <- cX_h  %>% filter(LD == 1) %>% filter(F_MISS <=0.01) %>% filter(HWE_P > (10**(-7)))
cXY_keep <- cXY_h  %>% filter(LD == 1) %>% filter(F_MISS <=0.01) %>% filter(HWE_P > (10**(-7)))

# write it out
vars.to.keep.all <- c(vars.to.keep$V1, cX_keep$SNP, cXY_keep$SNP)
write_tsv(data.frame(vars.to.keep.all), sprintf("%s/snp_filt_list_wX_v3.txt", DATA.DIR), col_names=FALSE)



