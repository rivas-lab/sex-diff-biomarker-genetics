# retrieve_sex_labels.R
# E Flynn
# 9/11/2017, updated 11/17/17 to downsample, 10/13/19 to include menopause labels
#
# Code to retrieve sex labels for IDs. 
# Retrieves sex labels and then filters to remove 151k individuals not included in the analysis. 

COVARIATE_MATRIX <- '/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
PHE_OUT_DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles'

# read in covariate matrix
cov_mat <- read.table(COVARIATE_MATRIX, header=TRUE)
cov_mat_sm <- cov_mat[,c("IID", "sex")]

# Filter out 151k sqc individuals
removal_file_one <- '/oak/stanford/groups/mrivas/ukbb/24983/sqc/ukb24983_remove.phe'
removal_file_two <- '/oak/stanford/groups/mrivas/ukbb/24983/sqc/w24983_20181016.csv'
remove1 <- read.table(removal_file_one, header=FALSE)
remove2 <- read.table(removal_file_two, header=FALSE) 
#ids.to.remove <- c(remove1[,1], unique(remove2[,1]))
#cov_mat_filt <- cov_mat_sm[!(cov_mat_sm$IID %in% ids.to.remove),]
cov_mat_filt <- cov_mat_sm

# zerosex
zeros <- cov_mat_filt[cov_mat_filt$sex==0,]
zeros <- zeros[zeros$IID > 0, ] ### remove negative IIDs



# onesex 
ones <- cov_mat_filt[cov_mat_filt$sex==1,]
ones <- ones[ones$IID > 0,] ### remove negative IIDs


# training/valid split
set.seed(1021)
num.f <- nrow(zeros) #181064
num.m <- nrow(ones) #156135

num.train.f <- floor(0.7*num.f)
num.train.m <- floor(0.7*num.m)

f_permuted_idx <- sample(1:num.f, num.f, replace=FALSE)
m_permuted_idx <- sample(1:num.m, num.m, replace=FALSE)

f_train_idx <- f_permuted_idx[1:num.train.f]
f_valid_idx <- f_permuted_idx[(num.train.f + 1):num.f]
rownames(zeros) <- NULL
zeros_train <- zeros[f_train_idx,]
zeros_valid <- zeros[f_valid_idx,]
rownames(ones) <- NULL
m_train_idx <- m_permuted_idx[1:num.train.m]
m_valid_idx <- m_permuted_idx[(num.train.m + 1):num.m]
ones_train <- ones[m_train_idx,]
ones_valid <- ones[m_valid_idx,]



# get pre vs post menopause
BASE.DIR <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas"
#meno_phe <- read.table(sprintf("%s/phe_extraction/menopause_phe_table.txt", BASE.DIR), header=TRUE)
#pre_meno <- meno_phe[meno_phe$meno.label=="pre",]
#post_meno <- meno_phe[meno_phe$meno.label=="post",]

# filter - all must be in the zeros!
#pre_meno2 <- pre_meno[pre_meno$IID %in% zeros$IID,]
#post_meno2 <- post_meno[post_meno$IID %in% zeros$IID,]


### DOWNSAMPLE
#num.pre <- nrow(pre_meno2)
#num.post <- nrow(post_meno2)


set.seed(1117)
keep.rows <- sample(1:num.f, num.m, replace=FALSE)
zeros.down <- zeros[keep.rows,]

#keep.rows2 <- sample(1:num.post, num.pre, replace=FALSE)
#keep.rows3 <- sample(1:num.m, num.pre, replace=FALSE)
#post.down <- post_meno2[keep.rows2,]
#ones.down <- ones[keep.rows3,]

# write it out
write.table(ones$IID, file=sprintf("%s/onesex.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)  	
write.table(zeros$IID, file=sprintf("%s/zerosex.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	

write.table(ones_train$IID, file=sprintf("%s/onesex_train.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)  	
write.table(zeros_train$IID, file=sprintf("%s/zerosex_train.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	

write.table(ones_valid$IID, file=sprintf("%s/onesex_valid.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)  	
write.table(zeros_valid$IID, file=sprintf("%s/zerosex_valid.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	


write.table(pre_meno2$IID, file=sprintf("%s/pre_meno.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	
write.table(post_meno2$IID, file=sprintf("%s/post_meno.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	

write.table(zeros.down$IID, file=sprintf("%s/zerosex_d.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)	
write.table(ones.down$IID, file=sprintf("%s/onesex_d.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(post.down$IID, file=sprintf("%s/post_meno_d.keep", PHE_OUT_DIR), col.names=FALSE, row.names=FALSE, quote=FALSE)		

