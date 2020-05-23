# extract_med_dat.R
# 11/11/2019
# Code for extracting medication data and making it tidy


require('tidyverse')
require('data.table')

med_visit_dat <- fread("/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/hormone_med/01_extract/med_extract.tsv", data.table=FALSE)

# combine columns together
# sadily the na.rm does not work even w the latest tidyr
united_dat <-
med_visit_dat %>% 
unite( "0", f.20003.0.0:f.20003.0.47, sep = ";", na.rm=TRUE)  %>%  
unite( "1", f.20003.1.0:f.20003.1.47, sep = ";", na.rm=TRUE) %>%
unite( "2", f.20003.2.0:f.20003.2.47, sep = ";", na.rm=TRUE) 

# get the visit/med setup 
long_united <- united_dat %>% pivot_longer(cols = v0:v2, names_to="visit", values_to="med_code")

# remove nas
long_united_narm <- long_united %>% mutate(med_code=gsub(";NA", "", med_code)) 
long_united_no_na <- long_united_narm %>% filter(med_code!="NA")

# separate rows
long_tidy <- long_united_no_na %>% separate_rows(med_code, sep=";") 

# write it out
write_tsv(long_tidy, "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/hormone_med/02_tidy/med_dat_visit_tidy.tsv")