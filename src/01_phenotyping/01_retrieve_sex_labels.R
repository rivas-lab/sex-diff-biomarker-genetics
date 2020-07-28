# 01_retrieve_sex_labels.R
# E Flynn
# shortened 05/2020
#
# Code to retrieve sex labels for IDs. 
# Retrieves sex labels and then filters to remove 151k individuals not included in the analysis. 

COVARIATE_MATRIX <- 'ukb24983_GWAS_covar.phe'
PHE_OUT_DIR <- 'phefiles'

# read in covariate matrix
cov_mat <- read.table(COVARIATE_MATRIX, header=TRUE)
cov_mat_sm <- cov_mat[,c("IID", "sex")]

# Filter out 151k sqc individuals
removal_file_one <- 'ukb24983_remove.phe'
removal_file_two <- 'w24983_20181016.csv'
remove1 <- read.table(removal_file_one, header=FALSE)
remove2 <- read.table(removal_file_two, header=FALSE) 

cov_mat_filt <- cov_mat_sm

# zerosex
zeros <- cov_mat_filt[cov_mat_filt$sex==0,]
zeros <- zeros[zeros$IID > 0, ] ### remove negative IIDs

# onesex 
ones <- cov_mat_filt[cov_mat_filt$sex==1,]
ones <- ones[ones$IID > 0,] ### remove negative IIDs

# write it out
write.table(ones$IID, file=sprintf("%s/onesex.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)  	
write.table(zeros$IID, file=sprintf("%s/zerosex.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	



