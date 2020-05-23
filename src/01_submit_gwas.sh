#!/bin/bash


# submit GWAS
GWAS_DIR=data/gwas_0522/
mkdir -p $GWAS_DIR

# ----- submit anthro traits ----- #
anthro_traits=(whr lfr afr tfr)
for trait in "${anthro_traits[@]}"; do
  sbatch src/03_gwas/run_sex_div_gwas.sh $trait zerosex phefiles/anthro/ $GWAS_DIR;
  sbatch src/03_gwas/run_sex_div_gwas.sh $trait onesex phefiles/anthro/ $GWAS_DIR;

    echo $trait;
done

# ------ submit biomarker traits ------ #
BIOMARKER_PHE_DIR=phefiles/biomarker_m_v2/
for phe_f in `ls $BIOMARKER_PHE_DIR`; do 
  trait="${phe_f%.*}";
  sbatch src/03_gwas/run_sex_div_gwas.sh $trait zerosex phefiles/biomarker_f_v2/ $GWAS_DIR;
  sbatch src/03_gwas/run_sex_div_gwas.sh $trait onesex phefiles/biomarker_m_v2/ $GWAS_DIR;
  echo $trait; 
done

BIOMARKER_PHE_DIR=phefiles/biomarker_m_v2/
for phe_f in `ls $BIOMARKER_PHE_DIR`; do 
  trait="${phe_f%.*}";
  echo $trait; 
done
