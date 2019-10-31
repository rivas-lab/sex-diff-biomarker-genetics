#
#
#

require('tidyverse')
require('dplyr')
require('reshape2')
options(stringsAsFactors=FALSE)


ss_phe <- read.table("../phe_extraction/sex_spec_pheno.txt", header=TRUE)
ss_phe_long <- melt(ss_phe, id.vars=c("IID", "age", "sex"))

ss_phe_long2 <- ss_phe_long %>% separate(variable, c(NA, "trait", "visit", NA))
ss_mat <- dcast(ss_phe_long2,  IID + visit ~ trait, value.var="value")

keep.rows <- apply(ss_mat[,3:ncol(ss_mat)], 1, function(x) any(!is.na(x)))
table(keep.rows)
ss_mat2 <- ss_mat[keep.rows,]


ss_mat3 <- filter(ss_mat2, IID >= 0)


ss_mat4 <- full_join(select(ss_phe, IID, sex, age), ss_mat3)

write.table(ss_mat4, file="../data/sex_spec_factor_mat.txt", row.names=FALSE, quote=FALSE, sep="\t")
