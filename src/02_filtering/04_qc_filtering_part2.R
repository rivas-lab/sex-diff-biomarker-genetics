require('tidyverse')
require('data.table')
DATA.DIR <- "data/"

# read in a file
options(stringsAsFactors=FALSE)



options(stringsAsFactors=FALSE)
vars.to.keep <- read.table(sprintf("%s/snp_filt_list.txt" , DATA.DIR), header=FALSE)
head(vars.to.keep)




# read in all the chr qc
cX_ld <- read.table(sprintf("%s/chr_qc/chrX_ld.prune.in", DATA.DIR), header=FALSE)
cXY_ld <- read.table(sprintf("%s/chr_qc/chrXY_ld.prune.in", DATA.DIR), header=FALSE)
cX_tab <- read.table(sprintf("%s/chr_qc/chrX_qc_table.txt", DATA.DIR), sep=" ", header=TRUE)
cXY_tab <- read.table(sprintf("%s/chr_qc/chrXY_qc_table.txt", DATA.DIR), sep=" ", header=TRUE)
cXY_hardy <- read_tsv(sprintf("%s/chr_qc/chrXY_hardy.hardy", DATA.DIR))
cX_hardy <- read_tsv(sprintf("%s/chr_qc/chrX_hardy.hardy.x", DATA.DIR))

intersect(cX_tab$SNP,vars.to.keep$V1 )
head(cX_tab)

table(cX_tab$keep==1)
table(cXY_tab$keep==1)

head(cX_tab)
head(cXY_tab)

filter(cX_tab, SNP %in% cX_ld$V1) %>% filter(LD !=1) %>% nrow()
filter(cXY_tab, SNP %in% cXY_ld$V1) %>% filter(LD !=1) %>% nrow()


filter(cX_tab, SNP %in% cX_ld$V1) %>% filter(LD ==1) %>% nrow()
filter(cXY_tab, SNP %in% cXY_ld$V1) %>% filter(LD==1) %>% nrow()


combinedX <- cX_tab %>% mutate(LD2=ifelse(SNP %in% cX_ld$V1,1,0))
combinedXY <- cXY_tab %>% mutate(LD2=ifelse(SNP %in% cXY_ld$V1,1,0))


combinedX2 <- combinedX %>% select(-LD, -keep) %>% rename(LD=LD2) 
combinedXY2 <- combinedXY %>% select(-LD, -keep) %>% rename(LD=LD2) 
head(combinedX2)
head(combinedXY2)

head(cX_hardy)
head(cXY_hardy)

cXY_h <- left_join(combinedXY2, cXY_hardy %>% select(ID, P), by=c("SNP"="ID")) %>% rename(HWE_P=P)
head(cXY_h)

cX_h <- left_join(combinedX2, cX_hardy %>% select(ID, P), by=c("SNP"="ID")) %>% rename(HWE_P=P)
head(cX_h)
nrow(cX_h)
nrow(cXY_h)

# this
cX_keep <- cX_h  %>% filter(LD == 1) %>% filter(F_MISS <=0.01) %>% filter(HWE_P > (10**(-7)))
cXY_keep <- cXY_h  %>% filter(LD == 1) %>% filter(F_MISS <=0.01) %>% filter(HWE_P > (10**(-7)))
#cX_keep <- cX_h  %>% filter(LD == 1) %>% filter(MAF >= 0.01) %>% filter(F_MISS <=0.1) %>% filter(HWE_P <= (10**(-7)))
#cXY_keep <- cXY_h  %>% filter(LD == 1) %>% filter(MAF >= 0.01) %>% filter(F_MISS <=0.1) %>% filter(HWE_P <= (10**(-7)))


#nrow(cX_keep)
#nrow(cXY_keep)

#head(cX_keep)

nrow(cX_keep)
nrow(cXY_keep)

head(cX_keep)

vars.to.keep.all <- c(vars.to.keep$V1, cX_keep$SNP, cXY_keep$SNP)


vars.to.keep.all <- c(vars.to.keep$V1, cX_keep$SNP, cXY_keep$SNP)
head(data.frame(vars.to.keep.all))

write_tsv(data.frame(vars.to.keep.all), sprintf("%s/snp_filt_list_wX_v3.txt", DATA.DIR), col_names=FALSE)

length(vars.to.keep.all)


