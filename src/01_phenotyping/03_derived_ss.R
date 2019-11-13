
require('tidyverse')
require('data.table')

# calculate the derived ss phenotypes
ss <- read_tsv("../data/hormone_med/02_tidy/ss_tidy.tsv")
age <- read_tsv("../data/hormone_med/02_tidy/age_tidy.tsv")

age_ss <- left_join(age, ss, by=c("f.eid", "visit")) 
# TODO - this removes some missing data  (e.g. no visit/no age) - yay, but we lose those 73 IDs

# add in sex
# TODO - actually extract this
COVARIATE_MATRIX <- '/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
cov_mat <- fread(COVARIATE_MATRIX, data.table=FALSE)
sex_lab <- cov_mat %>% select("IID", "sex") %>% rename("f.eid"="IID")

ss_age_sex <- full_join(sex_lab, age_ss)


# now derive
int_cols <- colnames(ss_age_sex)[5:ncol(ss_age_sex)]
colnames(ss_age_sex)[5:ncol(ss_age_sex)] <- sapply(int_cols, function(x) sprintf("X%s", x))
    
ss_mat <- ss_age_sex %>% mutate(
    preg=ifelse(X3140==1 | X3140==2, 1, X3140),
    pill=ifelse(X2804==-11, 1, ifelse(X2804 == -3 | X2804==-1,-9, 0)),
    hrt=ifelse(X3546==-11, 1, ifelse(X3546 == -3 | X3546==-1, -9, 0)))
    
# add in menopause data
ss_mat2 <- ss_mat %>% mutate(
    ooph= ifelse(X2834 == -3 | X2834==-5,-9, X2834),
    hyster= ifelse(X3591 == -3 | X3591==-5,-9, X3591),
    meno = ifelse(X2724 ==-3 | X2724==3,NA, ifelse(X2724 == 2 | X2724== 1, 1, X2724)), # hysterectomy / yes --> yes; not sure/prefer not --> -9
    hyster2= ifelse(X2724==2, 1, 0),
    meno.age=ifelse(X3581 == -3 | X3581==-1,-9, X3581))

ss_mat3 <- mutate(ss_mat2, 
                  years.post=ifelse(meno.age<0, NA, age - meno.age),# compute the years since menopause
                  surgical.meno=ifelse(ooph == 1 | hyster==1 | hyster2 == 1,1,ifelse(X2834 == -9 & X3591==-9, NA, 0))
                 ) 
    
# Period data
ss_mat4 <- ss_mat3 %>% mutate(
    period_today=ifelse(X3720 == -3 | X3720==-1,NA,X3720),
    irregular=ifelse(X3710 == -3 | X3710==-1,NA, ifelse(X3710==-6, 1, 0)), # -6 --> irregular
    cycle_length = ifelse(X3710 < 12, NA, ifelse(X3710 > 60, NA, X3710)) # filter > 60d
    ) %>% mutate(
    outside_cycle =ifelse(is.na(cycle_length), NA, ifelse(X3700 > 1.5*cycle_length, 1, 0))) %>% mutate(
        day_in_cycle = ifelse(X3700 <0 , NA, ifelse(X3700 > 60 | 
                                                    (!is.na(outside_cycle) & outside_cycle==1), NA, X3700)))
ss_mat5 <- ss_mat4 %>% mutate(
    normalized_day_in_cycle=ifelse(irregular | is.na(cycle_length) | is.na(day_in_cycle), NA,
                                   day_in_cycle/cycle_length*28)) %>% mutate(
    menstrual_phase=ifelse(period_today, "menstrual",
                            ifelse(is.na(normalized_day_in_cycle) | normalized_day_in_cycle > 15 & period_today == 1, NA, 
                                   ifelse(normalized_day_in_cycle < 10, "follicular",
                                        ifelse(normalized_day_in_cycle < 17, "peri-ov",
                                            ifelse(normalized_day_in_cycle < 26, "luteal", "late-luteal"))))))

# remove ID columns and write out
num_cols <- colnames(ss_mat5)[sapply(colnames(ss_mat5), function(x) substr(x, 0, 1)=="X")]
ss_mat6 <- ss_mat5 %>% select(- num_cols)
write_tsv(ss_mat6, "../../data/hormone_med/03_combine/ss_deriv.tsv")