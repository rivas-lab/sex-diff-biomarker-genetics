#!/bin/bash
#SBATCH --job-name=sex_div_gwas
#SBATCH --output=logs/sex_div_gwas_%A_%a.out
#SBATCH --error=logs/sex_div_gwas_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH -p rbaltman
#SBATCH --mem=10000


SEX=$1

PHE_ID=$2
if [ "$#" -ne 3 ]; then
  
  PHE_FILE='/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles/'${PHE_ID}'.phe'
else 
	PHE_FILE=$3
fi

CHR=$SLURM_ARRAY_TASK_ID

echo $PHE_FILE
echo $CHR

# for now some messy exports to set up dependencies, will fix later
export PATH=/oak/stanford/groups/mrivas/software/plink2/20170904/:$PATH
#export PATH=/home/erflynn/applications/plink2/20180223/:$PATH # location of latest plink version
export PATH=$PATH:$HOME/applications/htslib-1.5/installation/bin
ml py-scipystack # for pandas

./gwas_qt_sex_div_no_na.sh 24893 $SEX $PHE_ID $PHE_FILE $CHR



