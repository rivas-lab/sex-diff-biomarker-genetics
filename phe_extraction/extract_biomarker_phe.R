# extract_biomarker_phe.R
# E Flynn
# 5/20/2019 
#
# Code for extracting a biomarker phenotype from the matrix.
# 
# Usage: Rscript extract_biomarker_phe.R <phe_id>
# 
# note - requires 8MB memory

PHE_OUT_DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles'

args <- commandArgs(trailingOnly=TRUE)
phe_id <- args[1]

biomarker_phe <- read.delim("/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/phenotypes/raw/biomarkers_serum_full.phe")
print("read")
df_phe <- biomarker_phe[,c("f.eid", sprintf("f.%s.0.0", phe_id))]  # structured as phe, visit, info, we want visit 1
print("extracted")
df_phe2 <- cbind(df_phe[,"f.eid"], df_phe)
write.table(df_phe2, file=sprintf("%s/%s.phe", PHE_OUT_DIR, phe_id), row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")


