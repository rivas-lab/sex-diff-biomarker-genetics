#!/bin/bash
#
# 01_setup_traits.sh
# E Flynn
#
# This script cleans up the phe files for GWAS.


# 0. create the directories
mkdir -p phefiles/
mkdir -p phefiles/anthro
mkdir -p phefiles/biomarker_m
mkdir -p phefiles/biomarker_f

# 1. move phefiles to the directories
# ## FILL IN with your code
# these files have a .phe suffix and are formatted as "id id phenotype_value"

# compute + QC the phefiles
Rscript src/03_gwas/01_retrieve_sex_labels.R
Rscript src/03_gwas/02a_phe_qc_anthro.R 
Rscript src/03_gwas/02b_phe_qc_biomarker.R

