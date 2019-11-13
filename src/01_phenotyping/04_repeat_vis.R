# 04_repeat_vis.R
# 11/12/2019
# E Flynn
#
# Extract repeat visit data.
# We are interested in visits 0 and 1 where subjects have biomarker measurements.

require('tidyverse')
require('data.table')

# sex-specific data
ss_mat6 <- fread("../../data/hormone_med/03_combine/ss_deriv.tsv")
v1_ppl <- ss_mat6 %>% filter(visit==1)
v1_rows <- ss_mat6 %>% filter(f.eid %in% v1_ppl$f.eid & visit %in% c(0,1)) %>% arrange(f.eid) %>% mutate(f.eid=as.character(f.eid))

# biomarker data
biomarker_tidy <- fread("../../data/hormone_med/02_tidy/biomarker_tidy.txt", data.table=FALSE)
# remove the unimportant columns
phe_codes <- read_csv("../../data/ListPheCodes.csv")

biomarkers <- phe_codes %>% filter(category=="biomarker")
biomarker_cols <- intersect(biomarkers$trait, colnames(biomarker_tidy)) # missing creatinine in urine... 30510

biomarker_sm <- biomarker_tidy %>% select(IID, visit, biomarker_cols) %>% rename(f.eid=IID)


repeat_vis <- biomarker_sm %>% mutate(f.eid=as.character(f.eid)) %>% filter(visit==1)
repeat_vis %>% nrow()
length(intersect(repeat_vis$f.eid, v1_rows$f.eid )) # GREAT! consistency =)

biomarker_repeat <-biomarker_sm %>% filter(f.eid %in% v1_rows$f.eid)
ss_repeat <- v1_rows

# med data
med_dat <- fread("../../data/hormone_med/02_tidy/med_dat_visit_tidy.tsv", data.table=FALSE) %>%
mutate(f.eid=as.character(f.eid)) %>%
mutate(visit=as.numeric(str_replace(visit, "v", "")) )
med_repeat <- med_dat %>% filter(f.eid %in% v1_rows$f.eid) %>% filter(visit %in% c(0,1))

med_repeat %>% write_tsv("../../data/hormone_med/04_repeat_vis/med_repeat.tsv")
ss_repeat %>% write_tsv("../../data/hormone_med/04_repeat_vis/ss_repeat.tsv")
biomarker_repeat %>% write_tsv("../../data/hormone_med/04_repeat_vis/biomarker_repeat.tsv")