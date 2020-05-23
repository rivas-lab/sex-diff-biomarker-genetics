#!/bin/bash
# This script sets up the phe files to run GWAS


# create the directories
mkdir -p phefiles/
mkdir -p phefiles/anthro
mkdir -p phefiles/biomarker_m
mkdir -p phefiles/biomarker_f
mkdir -p phefiles/biomarker_m_v2
mkdir -p phefiles/biomarker_f_v2

# copy over the phefiles

# -- anthro -- #
cp /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/24611/phe/INI48.phe phefiles/anthro/
cp /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/24611/phe/INI49.phe phefiles/anthro/


# -- biomarker -- #
cp `/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/outputExtendedBMIreducedMaleWhiteBritish/phenotypes/residual/*.phe` phefiles/biomarker_m/

cp `/oak/stanford/groups/mrivas/projects/biomarkers/covariate_corrected/outputExtendedBMIreducedFemaleWhiteBritish/phenotypes/residual/*.phe` phefiles/biomarker_f/

# move the adjust statins to the actual traits
mv phefiles/biomarker_f/Apolipoprotein_B.adjust.statins.phe phefiles/biomarker_f/Apolipoprotein_B.phe
mv phefiles/biomarker_f/Cholesterol.adjust.statins.phe phefiles/biomarker_f/Cholesterol.phe
mv phefiles/biomarker_f/LDL_direct.adjust.statins.phe phefiles/biomarker_f/LDL_direct.phe
mv phefiles/biomarker_m/Apolipoprotein_B.adjust.statins.phe phefiles/biomarker_m/Apolipoprotein_B.phe
mv phefiles/biomarker_m/Cholesterol.adjust.statins.phe phefiles/biomarker_m/Cholesterol.phe
mv phefiles/biomarker_m/LDL_direct.adjust.statins.phe phefiles/biomarker_m/LDL_direct.phe

# compute + QC the phefiles
Rscript src/03_gwas/01_retrieve_sex_labels.R
Rscript src/03_gwas/02a_phe_qc_anthro.R 
Rscript src/03_gwas/02b_phe_qc_biomarker.R

