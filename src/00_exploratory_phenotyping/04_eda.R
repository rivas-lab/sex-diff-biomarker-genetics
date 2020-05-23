
require('tidyverse')
require('data.table')


# ---- read in and meno label ---- #
med_repeat <- read_tsv("../data/hormone_med/04_repeat_vis/med_repeat.tsv")
ss_repeat <- read_tsv("../data/hormone_med/04_repeat_vis/ss_repeat.tsv")
biomarker_repeat <- read_tsv("../data/hormone_med/04_repeat_vis/biomarker_repeat.tsv")



label_col <- function(meno, meno.age, age, sex, years.post, surgical.meno){
(if (sex==1 & !is.na(sex)){
    "male"
} else if (any(sapply(c(meno, age, sex), is.na))){
    "missing"
} else if (surgical.meno==1){
    "surgical_meno"
} else if (meno==0){
    ifelse(age > 60, "likely_meno", "pre")
} else {
    if (is.na(meno.age)) {
        "missing_age"
    } else{
           ifelse(meno.age <= 40, "premature",
           ifelse(years.post < 2 | is.na(years.post), "peri", "post")
          ) 
    }


})}
 


df <- ss_repeat
df2 <- df
df2$meno.label <- mapply(label_col, df$meno, df$meno.age, df$age, df$sex, df$years.post, df$surgical.meno)
meno_df <- df2 %>% select(f.eid, sex, age, visit, preg, pill, hrt, ooph, hyster, meno, hyster2, meno.age, years.post, surgical.meno, meno.label)


# ---- Nasa's code for setting up a matrix of covariates ---- #

library(dplyr)
anthro = data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/phenotypes/anthropometrics.phe", header=T, stringsAsFactors=F)
cov = data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/output/covariates/full.covariates.phe", header=T, stringsAsFactors=F)
alcohol <- data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/covariates/alcohol.phe", header=T, stringsAsFactors = F
)
smoking <- data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/covariates/smoking_status.phe", header=T, stringsAsFactors = F)
age <- data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/phenotypes/age_at_attendance.phe", header=T, stringsAsFactors = F)

statins <- data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/covariates/statins.phe", header=T, stringsAsFactors=F)
statins_followup <- data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/covariates/statins_followup.phe", header=T, stringsAsFactors=F)

marker.list <- read.table("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/names.txt", header=F, stringsAsFactors = F)[1:34,]

serum = data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/phenotypes/raw/biomarkers_serum_full.phe", header=T, stringsAsFactors=F)
urine = data.table::fread("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/phenotypes/raw/biomarkers_urine_full.phe", header=T, stringsAsFactors=F)

biomarkers <- inner_join(serum, urine)

template <- "f.%s.%d.0"

initial <- cov %>% select(FID, IID, sex, Array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                                                 PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20,
                           TDI) %>% inner_join(statins, by=c("IID" = "f.eid"))


initial.biomarkers <- marker.list$V2
names(initial.biomarkers) <- sprintf(template, marker.list$V1, 0)

# filter b/c some problems if I dont
initial.biomarkers <- initial.biomarkers[names(initial.biomarkers) %in% colnames(biomarkers)]


kept.initial <- biomarkers[,c("f.eid", names(initial.biomarkers))]
colnames(kept.initial) <- c("IID", initial.biomarkers)

initial <- inner_join(initial, kept.initial)
initial <- inner_join(initial, age %>% mutate(IID=f.eid, age=f.21003.0.0) %>% select(IID, age))
initial <- inner_join(initial, anthro %>% mutate(IID=f.eid,
                                                 weight=f.21002.0.0,
                                                 BMI=f.21002.0.0*10000/f.50.0.0**2,
                                                 WC=f.48.0.0,
                                                 HC=f.49.0.0,
                                                 WHR=WC/HC) %>% select(IID, weight, BMI, WHR, WC, HC))

initial <- inner_join(initial, alcohol %>% mutate(IID=f.eid, alcohol=ifelse(f.1618.0.0 == -6, 0.5,
                                                                            ifelse(f.1618.0.0 < 0, NA,
                                                                                   f.1618.0.0))) %>% select(IID, alcohol))

initial <- inner_join(initial, smoking %>%
                         mutate(IID=f.eid, smoking=ifelse(f.20116.0.0 == -3, NA, f.20116.0.0)) %>%
                         select(IID, smoking))

followup <- cov %>% select(FID, IID, sex, Array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                                                 PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20,
                           TDI) %>% inner_join(statins_followup, by=c("IID" = "f.eid"))

followup.biomarkers <- marker.list$V2
names(followup.biomarkers) <- sprintf(template, marker.list$V1, 1)

# filter b/c some problems if I dont
followup.biomarkers <- followup.biomarkers[names(followup.biomarkers) %in% colnames(biomarkers)]
print(followup.biomarkers)


kept.followup <- biomarkers[,c("f.eid", names(followup.biomarkers))]
colnames(kept.followup) <- c("IID", followup.biomarkers)

followup <- inner_join(followup, kept.followup)

followup <- inner_join(followup, age %>% mutate(IID=f.eid, age=f.21003.1.0) %>% select(IID, age))
followup <- inner_join(followup, anthro %>% mutate(IID=f.eid,
                                                 weight=f.21002.1.0,
                                                 BMI=f.21002.1.0*10000/f.50.1.0**2,
                                                 WC=f.48.1.0,
                                                 HC=f.49.1.0,
                                                 WHR=WC/HC) %>% select(IID, weight, BMI, WHR, WC, HC))

followup <- inner_join(followup, alcohol %>% mutate(IID=f.eid, alcohol=ifelse(f.1618.1.0 == -6, 0.5,
                                                                            ifelse(f.1618.1.0 < 0, NA,
                                                                                   f.1618.1.0))) %>% select(IID, alcohol))

followup <- inner_join(followup, smoking %>%
                         mutate(IID=f.eid, smoking=ifelse(f.20116.1.0 == -3, NA, f.20116.1.0)) %>%
                         select(IID, smoking))


combined <- inner_join(initial, followup, by=c("FID" = "FID", "IID" = "IID"), suffix=c(".initial", ".followup")) %>%
  inner_join(cov %>% select(FID, IID, sex, Array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                                                 PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20,
                           TDI))

combined$newStatins <- combined$AnyDrug.followup * !combined$AnyDrug.initial

combined$ageDiff <- combined$age.followup - combined$age.initial
combined$bmiDiff <- combined$BMI.followup - combined$BMI.initial
combined$whrDiff <- combined$WHR.followup - combined$BMI.initial

# ---- what nasa did before - examine how statins affect ----- #
combined.glms <- data.table::rbindlist(lapply(unlist(followup.biomarkers), function(trait) {
  trait.value.initial <- combined[,paste0(trait, ".initial")]
  trait.value.followup <- combined[,paste0(trait, ".followup")]
  tidied <- broom::tidy(glm(log(trait.value.followup) - log(trait.value.initial) ~ TDI + AnyDrug.initial + 
                            age.initial + ageDiff +
                BMI.initial + bmiDiff + age.initial*sex + ageDiff*sex + bmiDiff*sex +
                PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                newStatins + AnyDrug.initial + AnyDrug.followup, combined, family="gaussian"))
  tidied$trait <- trait
  tidied
}))
statin.effect <- combined.glms[combined.glms$term == "newStatins"]
statin.effect$q.value <- p.adjust(statin.effect$p.value, "fdr")
statin.effect[statin.effect$q.value < 0.01]


# ---- fit a glm for age/sex to examine covariates ---- #
age.sex.glms <- data.table::rbindlist(lapply(unlist(followup.biomarkers), function(trait) {
  trait.value.initial <- combined[,paste0(trait, ".initial")]
  trait.value.followup <- combined[,paste0(trait, ".followup")]
  tidied <- broom::tidy(glm(log(trait.value.initial) ~ TDI  + 
                            age.initial  +
                BMI.initial + sex + age.initial*sex + meno.initial +
                PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, combined, family="gaussian"))
  tidied$trait <- trait
  tidied
}))



age.sex.glms %>% filter(term %in% c( "age.initial","sex", "age.initial:sex")) %>%
mutate(q.value=p.adjust(p.value, "fdr")) %>% filter(q.value < 0.05) %>% write_tsv("../data/hormone_med/04_repeat_vis/age_sex.tsv")


# ---- reformat menopause data so that we can examine ---- #

meno_df2 <- meno_df %>% 
filter(meno.label %in% c("pre", "post"))  %>% 
filter(preg==0) %>% 
select(f.eid, age, visit, pill, hrt, meno.label) %>%
mutate(hrt=ifelse(is.na(hrt), 0, hrt), pill=ifelse(is.na(pill), 0, pill)) %>%
mutate(visit=ifelse(visit==0, "initial", "followup")) %>%
mutate(meno.label=ifelse(meno.label=="post", 1, 0)) 

m <- meno_df2 %>% pivot_wider(id_cols=f.eid, names_from=visit, names_prefix="meno.", values_from="meno.label") 
p <- meno_df2 %>% pivot_wider(id_cols=f.eid, names_from=visit, names_prefix="pill.", values_from="pill")
hrt <- meno_df2 %>% pivot_wider(id_cols=f.eid, names_from=visit, names_prefix="hrt.", values_from="hrt")

meno_df_wide <- inner_join(inner_join(m, p), hrt)
meno_df_wide <- meno_df_wide %>% 
mutate(pill.followup=ifelse(pill.followup==1, 1, 0),
       pill.initial=ifelse(pill.initial==1, 1, 0),
       hrt.initial=ifelse(hrt.initial==1, 1, 0),
       hrt.followup=ifelse(hrt.followup==1, 1, 0)) %>%
mutate(meno.change=ifelse(meno.initial==0 & meno.followup==1, 1, 0),
        pill.change=pill.followup-pill.initial,
       hrt.change=hrt.followup-hrt.initial
      )

meno_comb <- left_join(combined, meno_df_wide, by=c("IID"="f.eid"))
table(meno_df_wide$pill.change)
table(meno_df_wide$hrt.change)


# ---- examine menopause, pill, and hrt associations ---- #
meno.combined.glms <- data.table::rbindlist(lapply(unlist(followup.biomarkers), function(trait) {
  trait.value.initial <- meno_comb[,paste0(trait, ".initial")]
  trait.value.followup <- meno_comb[,paste0(trait, ".followup")]
  tidied <- broom::tidy(glm(log(trait.value.followup) - log(trait.value.initial) ~ 
                            TDI + age.initial + ageDiff +
                BMI.initial + bmiDiff + 
                PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                newStatins +  meno.change +
                             pill.change + hrt.change,
                            meno_comb, family="gaussian"))
  tidied$trait <- trait
  tidied
}))


meno.combined.glms %>%
mutate(q.value=p.adjust(p.value, "fdr")) %>% write_tsv("../data/hormone_med//04_repeat_vis/meno_fit.tsv")

term_est <- meno.combined.glms %>%
mutate(q.value=p.adjust(p.value, "fdr")) %>% 
filter(q.value < 0.05) %>% filter(term %in% colnames(meno_df_wide))



# ---- now let's look at hormone drugs!! ---- #
hormone_drugs <- read.delim("../data/hormone_drugs_todownload_v2.txt")
hormone_drugs_long <- hormone_drugs %>% 
select(drug, read_codes, hormone) %>% 
separate_rows(read_codes, sep=";") %>%
mutate(read_codes=str_trim(read_codes)) %>%
filter(read_codes!="")

read_codes <- hormone_drugs_long$read_codes
read_codes2 <- unique(read_codes)

med_repeat2 <- med_repeat %>% mutate(med_code=as.character(med_code))

rep_hormone_drug <- med_repeat2 %>% filter(med_code %in% read_codes2)
rep_hormone_drug <- left_join(rep_hormone_drug, hormone_drugs_long, by=c("med_code"="read_codes"))


# some of these increase hormones, others decrease
hormone_drugs %>% filter(hormone=="testosterone")
t_inc <- hormone_drugs %>% filter(drug %in% c("testosterone", "methyltestosterone")) %>% 
    separate_rows(read_codes, sep=";")
t_dec <- hormone_drugs %>% filter(drug=="spironolactone") %>% separate_rows(read_codes, sep=";")
hormone_drugs %>% filter(hormone=="estrogen")
e_inc <- hormone_drugs  %>% filter(drug %in% c("estradiol", "ethinyestradiol", "estriol")) %>% 
 separate_rows(read_codes, sep=";")
p_drugs <- hormone_drugs %>% filter(hormone=="progesterone") %>% separate_rows(read_codes, sep=";")

# ---- ID common meds ---- #
med_counts <- med_repeat %>% mutate(med_code=as.character(med_code)) %>% 
group_by(med_code) %>% 
count() %>% arrange(desc(n))  

coding <- read_tsv("/oak/stanford/groups/rbaltman/ukbiobank/metadata/coding4.tsv") %>% 
mutate(coding=as.character(coding))

med_counts_lab <- left_join(med_counts, coding, by=c("med_code"="coding")) %>% rename(drug=meaning)

common_drugs <- med_counts_lab %>% filter(n > 500)
# remove free text entry and multivitamins and statins
common_drugs2 <- common_drugs %>% filter(!med_code %in% c("99999", "1140852976")) %>% filter(! drug %in% c("simvastatin", "atorvastatin"))


common_repeats <- med_repeat %>% filter(med_code %in% common_drugs$med_code)
common_repeats$taken <- 1



# ---- functions for running the med analysis ----- #

# prepare the covariate data frame
prep_med_input <- function(med_dat, my_code, comb_tab){
    
    med <- med_dat %>% filter(med_code %in% c(my_code)) %>%
    mutate(med_code=1) %>%
    mutate(visit=ifelse(visit==0, "initial", "followup")) %>% 
    unique() %>%
    pivot_wider(id_cols=f.eid, names_from=visit, values_from=med_code, names_prefix="med.") 
head(med)
 med2 <- med %>% mutate(
    med.initial=ifelse(is.na(med.initial), 0, 1),
    med.followup=ifelse(is.na(med.followup), 0, 1)) %>%
mutate(
    med.start=ifelse(med.initial==0 & med.followup==1, 0, 1),
    med.stop=ifelse(med.initial==1 & med.followup==0, 0, 1))
    
left_join(combined, full_join(med2, comb_tab), by=c("IID"="f.eid")) %>%
    mutate(med.initial=ifelse(is.na(med.initial), 0, med.initial),
      med.followup=ifelse(is.na(med.followup), 0, med.followup),
      med.start=ifelse(is.na(med.start), 0, med.start),
          med.stop=ifelse(is.na(med.stop), 0, med.stop))
}


# function for fitting a given drug
get_drug_fit <- function(read_code, med_df, covar_df, drug_name=NULL){
    comb_med <- prep_med_input(med_df, read_code, covar_df)
    med_dat_other <- data.table::rbindlist(lapply(unlist(followup.biomarkers), function(trait) {
      trait.value.initial <- comb_med[,paste0(trait, ".initial")]
      trait.value.followup <- comb_med[,paste0(trait, ".followup")]
      tidied <- broom::tidy(glm(log(trait.value.followup) - log(trait.value.initial) ~ 
                                TDI + age.initial + ageDiff + age.initial*sex + ageDiff*sex + bmiDiff*sex +
                    BMI.initial + bmiDiff + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                    PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                    newStatins +  meno.change +
                                 pill.change + hrt.change+med.start + med.stop, 
                                comb_med, family="gaussian"))
      tidied$trait <- trait
      tidied
    })) 
    med_fit <- med_dat_other %>% filter(str_detect(term, "med")) %>%
    mutate(q.value=p.adjust(p.value, "fdr")) %>%  arrange(q.value) %>% filter(q.value < 0.05)
    if (is.null(drug_name)){drug_name=read_code[[1]]}
    med_fit$med <- rep(drug_name, nrow(med_fit))
    med_fit
}


# ---- fit all of the hormone drugs ---- #
t_inc_fit <- get_drug_fit(t_inc$read_codes, med_repeat, meno_df_wide, "inc_t_drugs")

t_dec_fit <- get_drug_fit(t_dec$read_codes, med_repeat, meno_df_wide, "dec_t_drugs")
e_inc_fit <- get_drug_fit(e_inc$read_codes, med_repeat, meno_df_wide, "inc_e_drugs")
p_drug_fit <- get_drug_fit(p_drugs$read_codes, med_repeat, meno_df_wide, "p_drugs")


hormone_drugs <- do.call(rbind, list(t_inc_fit, t_dec_fit, e_inc_fit, p_drug_fit))
hormone_drugs %>% write_tsv("../data/hormone_med/04_repeat_vis/hormone_drug_fit.tsv")


# ---- fit all of the common drugs ---- #
res <- do.call(rbind, lapply(unique(common_drugs2$med_code), function(drug) 
    get_drug_fit(drug, common_repeats %>% select(-taken), meno_df_wide)))

left_join(res, common_drugs %>% select(med_code, drug), by=c("med"="med_code")) %>% arrange(q.value) %>% 
write_tsv("../data/hormone_med/04_repeat_vis/drug_s_s_fit.tsv")


