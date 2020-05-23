require('tidyverse')
require('data.table')
PHE.DIR <- "../../../phefiles"
ANTHRO.DIR <- "../../../phefiles/anthro"


COVARIATE_MATRIX <- '/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
# read in covariate matrix
cov_mat <- fread(COVARIATE_MATRIX, data.table=FALSE) %>% select(IID, sex)

# function for pre-processing an anthro trait
#  - removes missing data (-9)
#  - removes those outside 6sd (separate for M+F)
processAnthroTrait <- function(phe){
    phe1 <- inner_join(phe, cov_mat, by=c("IID"))
    trait_name <- colnames(phe)[3]
    colnames(phe1)[3] <- "trait"
    phe2 <- phe1 %>% filter(trait>0)
    print(sprintf("Removed %s indiviuals who were missing info.", nrow(phe)-nrow(phe2)))
    phe2.1 <- phe2 %>% 
        group_by(sex) %>% 
        mutate(mu_s=mean(trait), s_s=sd(trait)) %>% 
        mutate(low_c=mu_s-(6*s_s),up_c=mu_s+(6*s_s))
    phe3 <- phe2.1 %>% ungroup() %>% 
        filter(trait > low_c & trait < up_c) %>% 
        select(FID, IID, trait)
    colnames(phe3)[3] <- trait_name
    print(sprintf("Removed %s indiviuals who were outside 6sd.", nrow(phe2)-nrow(phe3)))

    return(phe3)
}


# --- WHR --- #
wc <- read.table(sprintf("%s/INI48.phe", ANTRHO.DIR))
hc <- read.table(sprintf("%s/INI49.phe", ANTHRO.DIR))
colnames(wc) <- c("FID", "IID", "wc")
colnames(hc) <- c("FID", "IID", "hc")

# filter
hc_f <- processAnthroTrait(hc) %>% as_tibble()
wc_f <- processAnthroTrait(wc) %>% as_tibble()

# calculate
whr_tab <- inner_join(wc_f, hc_f %>% select(-FID), by=c("IID")) %>% mutate(whr=wc/hc)

write_tsv(whr_tab %>% select(-wc, -hc), sprintf("%s/whr.phe", ANTHRO.DIR))

# --- fat ratios --- #
total <- read.table(sprintf("%s/INI23100.phe", ANTHRO.DIR))
leg_l <- read.table(sprintf("%s/INI23116.phe", ANTHRO.DIR))
leg_r <- read.table(sprintf("%s/INI23112.phe", ANTHRO.DIR))
arm_l <- read.table(sprintf("%s/INI23124.phe", ANTHRO.DIR))
arm_r <- read.table(sprintf("%s/INI23120.phe", ANTHRO.DIR))
trunk <- read.table(sprintf("%s/INI23128.phe", ANTHRO.DIR))

colnames(leg_l) <- c("FID", "IID", "leg_l")
colnames(leg_r) <- c("FID", "IID", "leg_r")
colnames(arm_l) <- c("FID", "IID", "arm_l")
colnames(arm_r) <- c("FID", "IID", "arm_r")
colnames(trunk) <- c("FID", "IID", "trunk")
colnames(total) <- c("FID", "IID", "total")

# filter 
leg_l_f <- processAnthroTrait(leg_l) %>% as_tibble()
leg_r_f <- processAnthroTrait(leg_r) %>% as_tibble()
arm_l_f <- processAnthroTrait(arm_l) %>% as_tibble()
arm_r_f <- processAnthroTrait(arm_r) %>% as_tibble()
trunk_f <- processAnthroTrait(trunk) %>% as_tibble()
total_f <- processAnthroTrait(total) %>% as_tibble()

# calculate
legs <- inner_join(leg_l_f, leg_r_f %>% select(-FID), by="IID") %>% mutate(legs=leg_l+leg_r)
arms <- inner_join(arm_l_f, arm_r_f %>% select(-FID), by="IID") %>% mutate(arms=arm_l+arm_r)
legs2 <- inner_join(total_f, legs %>% select(-FID), by="IID") %>% mutate(lfr=legs/total)
arms2 <- inner_join(total_f, arms %>% select(-FID), by="IID") %>% mutate(afr=arms/total)
trunk2 <- inner_join(total_f, trunk_f %>% select(-FID), by="IID") %>% mutate(tfr=trunk/total)

# write it out
legs2 %>% select(FID, IID, lfr) %>% write_tsv(sprintf("%s/lfr.phe", ANTHRO.DIR))
arms2 %>% select(FID, IID, afr) %>% write_tsv(sprintf("%s/afr.phe", ANTHRO.DIR))
trunk2 %>% select(FID, IID, tfr) %>% write_tsv(sprintf("%s/tfr.phe", ANTHRO.DIR))