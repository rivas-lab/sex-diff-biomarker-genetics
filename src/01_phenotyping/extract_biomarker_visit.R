# extract_biomarker_visit.R
# E Flynn
# 6/13/2019
#
# Code for extracting biomarker data, and labeling with visit.
# The end result is a matrix that is IID, visit, biomarker_1...biomarker_n

require('tidyverse')
require('data.table')
require('reshape2')
options(stringsAsFactors=FALSE)


biomarker_phe <- fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/phenotypes/raw/biomarkers_serum_full.phe",
                      data.table=FALSE)
biomarker_phe <- rename(biomarker_phe, IID=f.eid)

# melt and separate the visit variables
bio_long <- melt(biomarker_phe, id.vars=c("IID"))
bio_long2 <- bio_long %>% separate(variable, c(NA, "trait", "visit", NA))

# put back into a matrix
bio_mat <- dcast(bio_long2,  IID + visit ~ trait, value.var="value")

# remove empty data 
keep.rows <- apply(bio_mat[,3:ncol(bio_mat)], 1, function(x) any(!is.na(x)))
table(keep.rows)
bio_mat2 <- bio_mat[keep.rows,]
table(bio_mat2$IID >=0)
bio_mat3 <- filter(bio_mat2, IID >= 0)

# write out results
write.table(bio_mat3, file="../data/biomarker_visit_mat.txt", row.names=FALSE, quote=FALSE)