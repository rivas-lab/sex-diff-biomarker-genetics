#!/bin/bash
#SBATCH --job-name=SS_GWAS
#SBATCH --output=logs/ss_gwas.%A_%a.out
#SBATCH --error=logs/ss_gwas.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=51200
#SBATCH --time=1:00:00

trait=$1
sex=$2
phe_dir=$3
out_dir=$4

ml load biology; ml load plink/1.90b5.3; ml load plink/2.0a2

python src/03_gwas/gwas.py --run-array --run-now --mem=51200 --cores=8 --pheno ${phe_dir}/${trait}.phe --out ${out_dir} --log-dir logs --sex-div True --keep-sex ${sex} --keep-sex-file phefiles/${sex}.keep --include-x True
