require('tidyverse')
PHE.DIR <- "../../../phefiles"
BIOMARKER_F_PATH <- "../../../phefiles/biomarker_f"
BIOMARKER_M_PATH <- "../../../phefiles/biomarker_m"

# get a list of traits
biomarker_files <- list.files(path=BIOMARKER_F_PATH, pattern="*.phe")
traits <- unique(sapply(biomarker_files, function(x) strsplit(x, "\\.")[[1]][[1]]))

traits.to.rem <- c("Fasting_glucose", "Oestradiol", "Rheumatoid_factor", "Microalbumin_in_urine")
list.traits <- data.frame(setdiff(traits, traits.to.rem), stringsAsFactors=FALSE)
colnames(list.traits) <- ""
write_tsv(list.traits, sprintf("%s/list_traits.txt", PHE.DIR))
                        
# function to remove NAs and rewrite it out in a new directory
remNaBio <- function(trait) {
    print(trait)
    df <- read_tsv(sprintf("%s/%s.phe", BIOMARKER_F_PATH, trait))
    df2 <- df[!is.na(df[,3]),] 
    df2 %>% write_tsv(sprintf("%s_v2/%s.phe", BIOMARKER_F_PATH, trait))
    print(nrow(df2)-nrow(df))

    rm(df)
    rm(df2)
    df <- read_tsv(sprintf("%s/%s.phe", BIOMARKER_M_PATH, trait))
    df2 <- df[!is.na(df[,3]),] 
    df2 %>% write_tsv(sprintf("%s_v2/%s.phe", BIOMARKER_M_PATH, trait))
        print(nrow(df2)-nrow(df))
}

# apply to list of traits
lapply(list.traits[,1], remNaBio)