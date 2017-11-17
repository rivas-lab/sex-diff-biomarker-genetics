

# exploratory analysis of xchr snps

DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"
QC.DIR <- sprintf('%s/chr_qc/', DATA.FOLDER)

source('model_utils.R')
require('qqman')
dat.f <- read.table("ukb24893_v2_cX.50.zerosex.PHENO1.glm.linear")
dat.m <- read.table("ukb24893_v2_cX.50.onesex.PHENO1.glm.linear")
col.labels <- c("CHROM", "POS", "ID", "REF", "ALT1", "TEST", "OBS_CT", 
        "BETA", "SE", "T_STAT", "P")
colnames(dat.f) <- col.labels
colnames(dat.m) <- col.labels
filt.dat <- filtUkbDat(dat.f, dat.m)

all.dat.f <- filt.dat$`1`
all.dat.m <- filt.dat$`2`
colnames(all.dat.f)[1:3] <- c("CHR", "BP", "SNP")
colnames(all.dat.m)[1:3] <- c("CHR", "BP", "SNP")
all.dat.f$CHR <- 23
all.dat.m$CHR <- 23
manhattan(all.dat.f, main="Height (women)")
manhattan(all.dat.m, main="Height (men)")


# se filt
nrow(all.dat.f) # 17521
se.filt <- filterSE(all.dat.f, all.dat.m, 'quant', 0.2) # 17930

xchr.qc <- read.table(sprintf("%s/chrX_qc_table.txt", QC.DIR))
vars.to.keep2 <- xchr.qc[xchr.qc$keep==1,]
# variant filter
#dim(se.filt$`1`[se.filt$`1`$SNP %in% vars.to.keep,])  ### 9396 
dim(se.filt$`1`[se.filt$`1`$SNP %in% vars.to.keep2,]) ### 9299
final.f <- se.filt$`1`[se.filt$`1`$SNP %in% vars.to.keep2,]
final.m <- se.filt$`2`[se.filt$`2`$SNP %in% vars.to.keep2,]
manhattan(final.f, main="Height (women)")
manhattan(final.m, main="Height (men)")
