# extractPhe.R
# E Flynn
# Last Updated: 8/11/2019

require('tidyverse')
require('data.table')
options(stringsAsFactors=FALSE)

PHENO.PATH <- "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata"
BASE.DIR <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas"

extractCode <- function(x) {
  strsplit(as.character(x), ".", fixed=TRUE)[[1]][[2]]
}

# --- pregnancy, age, and gender data --- #
UKB.F1 <- sprintf("%s/9796/21732/download/ukb21732.tab", PHENO.PATH)
metadata_info <- fread(UKB.F1, data.table=FALSE, nrows=1)

codes <- map(colnames(metadata_info), extractCode) %>% unlist()
metadata_codes <- which(codes %in% c("3140", "21003", "31", "22001")) # codes for pregnancy, age, gender, genetic sex respectively
metadata_str <- paste(c(1,metadata_codes), collapse=",")
system(sprintf("cut -f %s %s > %s/phe_extraction/baseline_char.txt", metadata_str, UKB.F1, BASE.DIR))
metadata <- fread(sprintf("%s/phe_extraction/baseline_char.txt", BASE.DIR), data.table=FALSE)
#metadata <- fread(UKB.F1, data.table=FALSE, select=c(1,metadata_codes))

# --- Sex-specific data --- #
UKB.F2 <- sprintf("%s/2000269/21730/download/ukb21730.tab", PHENO.PATH)
ssfile_cols <- fread(UKB.F2, data.table=FALSE, nrows=1)
codes_list <- map(colnames(ssfile_cols), extractCode) %>% unlist()

phe_codes <- read.csv(sprintf("%s/phe_extraction/ListPheCodes.csv", BASE.DIR)) # TODO - how did I get this list?
phe_codes$X <- NULL

sex_spec <- filter(phe_codes, category == "sex specific")
list_traits <- sex_spec$trait
ss_codes <- which(codes_list %in% list_traits)
system(sprintf("cut -f %s %s > %s/phe_extraction/ss_cols.txt", paste(c(1, ss_codes), collapse=","), UKB.F2, BASE.DIR))
ss_dat <- fread(sprintf("%s/phe_extraction/ss_cols.txt", BASE.DIR), data.table=FALSE)

# --- put the two together --- #


# --- filter - we only care abt visit #0 and visit #1 --- #

# ---- WRITE IT OUT