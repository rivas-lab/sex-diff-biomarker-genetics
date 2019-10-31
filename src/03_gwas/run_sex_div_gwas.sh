#!/bin/bash
#SBATCH --job-name=SS_GWAS
#SBATCH --output=logs/ss_gwas.%A_%a.out
#SBATCH --error=logs/ss_gwas.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=51200
#SBATCH --time=12:00:00
#SBATCH -p rbaltman,owners

trait=$1
sex=$2
phe_dir=$3
out_dir=$4

ml load biology; ml load plink; ml load plink2

BASE_DIR=/scratch/PI/mrivas/users/erflynn/sex_div_gwas


python gwas.py --run-array --run-now --mem=51200 --cores=8 --pheno ${phe_dir}/${trait}.phe --out ${out_dir} --log-dir ${BASE_DIR}/logs --sex-div True --keep-sex ${sex} --keep-sex-file ${BASE_DIR}/phefiles/${sex}.keep --include-x True
