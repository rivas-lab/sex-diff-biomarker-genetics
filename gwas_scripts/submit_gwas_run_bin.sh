#!/bin/bash
#SBATCH --job-name=sex_div_gwas
#SBATCH --output=sex_div_gwas_%j.out
#SBATCH --error=sex_div_gwas_%j.err
#SBATCH --time=2:00:00
#SBATCH -p normal
#SBATCH --mem=10000


SEX=$1

PHE_ID=$2
PHE_FILE='/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles/'${PHE_ID}'.phe'

CHR=$SLURM_ARRAY_TASK_ID

# for now some messy exports to set up dependencies, will fix later
export PATH=/oak/stanford/groups/mrivas/software/plink2/20170904/:$PATH
export PATH=$PATH:$HOME/applications/htslib-1.5/installation/bin
source activate root # loads my root env for pandas

./gwas_bin_sex_div_no_na.sh 24893 $SEX $PHE_ID $PHE_FILE $CHR 
