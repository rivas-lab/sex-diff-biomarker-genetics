require('tidyverse')
require('data.table')
ukb24611 <- fread("../../data/hormone_med/01_extract/ukb24611_age.tsv")

long_age <- ukb24611 %>%
 filter(f.eid > 0) %>%
 rename("0"=f.21003.0.0, "1"=f.21003.1.0, "2"=f.21003.2.0) %>%
 select(f.eid, "0", "1", "2") %>%
 pivot_longer(cols=c("0", "1", "2"), names_to="visit", values_to="age") %>% 
 filter(!is.na(age)) 

long_age %>% write_tsv("../../data/hormone_med/02_tidy/age_tidy.tsv")