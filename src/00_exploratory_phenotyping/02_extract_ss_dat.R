# extract_ss_dat.R
# 11/11/2019
# Code for extracting/tidying ss dat

require('tidyverse')
require('data.table')

ss_dat <- fread("../../data/hormone_med/01_extract/ss_dat.tsv", data.table=FALSE)

preg_dat <- fread("../../data/hormone_med/01_extract/preg_dat.tsv", data.table=FALSE)

ss_phe <- full_join(preg_dat, ss_dat, by="f.eid")

ss_phe_long <- melt(ss_phe, id.vars=c("f.eid"))

ss_phe_long2 <- ss_phe_long %>% separate(variable, c(NA, "trait", "visit", NA))
ss_mat <- dcast(ss_phe_long2,  f.eid + visit ~ trait, value.var="value")

write_tsv(ss_mat, "../../data/hormone_med/02_tidy/ss_tidy.tsv")
