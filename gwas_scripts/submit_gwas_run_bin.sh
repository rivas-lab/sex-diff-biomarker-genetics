#!/bin/bash
#SBATCH --job-name=sex_div_gwas
#SBATCH --output=logs/sex_div_gwas_%A_%a.out
#SBATCH --error=logs/sex_div_gwas_%A_%a.err
#SBATCH --time=3:00:00
#SBATCH -p rbaltman
#SBATCH --mem=20000


SEX=$1

PHE_ID=$2
PHE_FILE='/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles/'${PHE_ID}'.phe'

CHR=$SLURM_ARRAY_TASK_ID

# for now some messy exports to set up dependencies, will fix later
#export PATH=/oak/stanford/groups/mrivas/software/plink2/20170904/:$PATH
export PATH=/home/erflynn/applications/plink2/20180223/:$PATH # location of latest plink version
export PATH=$PATH:$HOME/applications/htslib-1.5/installation/bin
ml py-scipystack

./gwas_bin_sex_div_no_na.sh 24893 $SEX $PHE_ID $PHE_FILE $CHR 
