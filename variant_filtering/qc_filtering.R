# qc_filtering.R
# E Flynn
# 11/9/17

# Code for running initial QC filtering for X, Y, XY (PAR), MT chromosomes.
# Filters by LD, missingness <= 0.1, MAF >= 0.01. We do not have HWE data for X chromosome - this is tricky to get.


DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"
QC.DIR <- sprintf('%s/chr_qc/', DATA.FOLDER)
LD.DIR <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/gwas_scripts/ld_dat/"

generateQCFile <- function(chr){
	       # updated to more stringent LD filters to match other data
	ld <- read.table(sprintf("%s/ld_out%s.txt.prune.in", LD.DIR, chr), colClasses='character') # list of variants that pass  	
	#ld <- read.table(sprintf("%s/chr%s_ld.prune.in", QC.DIR, chr), colClasses='character') # list of variants that pass
	colnames(ld) <- c("SNP")

	f <- read.table(sprintf("%s/chr%s.f.afreq", QC.DIR, chr))
	colnames(f) <- c("CHR", "SNP", "A1", "A2", "MAF", "NOBS") # number of allele observations

	m <- read.table(sprintf("%s/chr%s.m.vmiss", QC.DIR, chr))
	colnames(m) <- c("CHR", "SNP", "N_MISS", "N_GENO", "F_MISS") # missing call rate

	combined <- merge(f[,c("CHR", "SNP", "A1", "A2", "MAF")], m[,c("SNP", "F_MISS")], by="SNP")
	combined$LD <- sapply(combined$SNP, function(x) ifelse(as.character(x) %in% ld$SNP, 1, 0))

	print(nrow(combined))
	#nrow(combined[combined$LD == 1 & combined$MAF >= 0.01,]) # 10297
	print(nrow(combined[combined$LD == 1 & combined$MAF >= 0.01 & combined$F_MISS <= 0.1,])) # 9849
	#vars.to.keep <- combined[combined$LD == 1 & combined$MAF >= 0.01,]$SNP
	vars.to.keep2 <- combined[combined$LD == 1 & combined$MAF >= 0.01 & combined$F_MISS <= 0.1,]$SNP

	combined$keep <- sapply(combined$SNP, function(x) ifelse(as.character(x) %in% vars.to.keep2, 1, 0))
	write.table(combined, sprintf("%s/chr%s_qc_table.txt", QC.DIR, chr), row.names=FALSE, quote=FALSE)
	return(combined)	
}

chr.tables <- lapply(c("XY", "X", "Y", "MT"), generateQCFile)
full.tab <- do.call(rbind, chr.tables)
write.table(full.tab, file="alt_chr_qc_table.txt", row.names=FALSE, quote=FALSE)


