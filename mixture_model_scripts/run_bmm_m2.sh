#!/bin/bash
#SBATCH --job-name=run_m2
#SBATCH --output=logs/run_m2_%a_%A.out
#SBATCH --error=logs/run_m2_%a_%A.err
#SBATCH --time=30:00:00
#SBATCH --partition=rbaltman
#SBATCH --array=1-2
#SBATCH --mem=50000

module load R/3.3.0

trait=$1
trait_type=$2
Rscript run_bmm_m2.R ${SLURM_ARRAY_TASK_ID} ${trait} ${trait_type}
