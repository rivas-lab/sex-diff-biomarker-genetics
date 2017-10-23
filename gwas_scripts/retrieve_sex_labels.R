# retrieve_sex_labels.R
# E Flynn
# 9/11/2017
#
# Code to retrieve sex labels for IDs. 
# Retrieves sex labels and then filters to remove 151k individuals not included in the analysis. 

COVARIATE_MATRIX <- '/scratch/PI/mrivas/ukbb/24983/phe/ukb24983_GWAS_covar.phe'
PHE_OUT_DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles'

# read in covariate matrix
cov_mat <- read.table(COVARIATE_MATRIX, header=TRUE)
cov_mat_sm <- cov_mat[,c("IID", "sex")]

# Filter out 151k sqc individuals
removal_file_one <- '/scratch/PI/mrivas/ukbb/24983/phe/ukb24983_remove.phe'
removal_file_two <- '/scratch/PI/mrivas/ukbb/24983/phe/w2498_20170726.phe'
remove1 <- read.table(removal_file_one, header=FALSE)
remove2 <- read.table(removal_file_two, header=FALSE) # this file is messy b/c of meta characters - only has three people in it tho
ids.to.remove <- c(remove1[,1], unique(remove2[,1]))
cov_mat_filt <- cov_mat_sm[!(cov_mat_sm$IID %in% ids.to.remove),]

# zerosex
zeros <- cov_mat_filt[cov_mat_filt$sex==0,]
zeros <- zeros[zeros$IID > 0, ] ### remove negative IIDs
write.table(zeros$IID, file=sprintf("%s/zerosex.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	

# onesex 
ones <- cov_mat_filt[cov_mat_filt$sex==1,]
ones <- ones[ones$IID > 0,] ### remove negative IIDs
write.table(ones$IID, file=sprintf("%s/onesex.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	



