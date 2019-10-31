# extract_hormonal_phe.R
# E Flynn
# 6/13/2019
#
# Code for extracting hormonal phenotypes.
# Note: code for extracting pregnancy phenotype is in BASH and is commented out. This must be run prior to starting

require('tidyverse')
require('data.table')
options(stringsAsFactors=FALSE)

phe_codes <- read.csv("../phe_extraction/ListPheCodes.csv")
phe_codes$X <- NULL

sex_spec <- filter(phe_codes, category == "sex specific")



## PREGNANCY DATA EXTRACTION ##
# head -1 /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/21732/download/ukb21732.tab > tmp_colnames.txt

# > preg_cols <- read.delim("tmp_colnames.txt")
# > preg_cols2 <- sapply(colnames(preg_cols), function(x) strsplit(as.character(x), ".", fixed=TRUE)[[1]][[2]]=="3140")
# > preg_cols2[[1]] <- TRUE
# > which(preg_cols2==TRUE) 

# cut -f 1,753,754,755 /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/21732/download/ukb21732.tab > ../phe_extraction/preg_data.txt
preg_dat <- fread("../phe_extraction/preg_data.txt")


tab_file <- fread('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2000269/21730/download/ukb21730.tab', data.table=FALSE)

list_cols <- colnames(tab_file)
cols_keep <- sapply(list_cols, function(x) strsplit(x, ".", fixed=TRUE)[[1]][[2]] %in% list_traits)


cols_keep[1] <- TRUE
dat <- tab_file[,cols_keep]
trait_counts <- table(sapply(colnames(dat), function(x) strsplit(x, ".", fixed=TRUE)[[1]][[2]] )) # three visits for each of nine traits


dat2 <- full_join(preg_dat, dat)

COVARIATE_MATRIX <- '/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
cov_mat <- read.table(COVARIATE_MATRIX, header=TRUE, stringsAsFactors=FALSE)
covar_data <- cov_mat[,c("IID","age", "sex")]
dat3 <- full_join(covar_data, dat2, c("IID"="f.eid"))

write.table(dat3, file="../phe_extraction/sex_spec_pheno.txt", row.names=FALSE, quote=FALSE, sep="\t")