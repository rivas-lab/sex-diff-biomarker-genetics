#!/bin/bash
#SBATCH --job-name=run_m2
#SBATCH --output=logs/run_m2_%a_%A.out
#SBATCH --error=logs/run_m2_%a_%A.err
#SBATCH --time=48:00:00
#SBATCH --partition=rbaltman
#SBATCH --array=2
#SBATCH --mem=50000

ml r-rstan

trait=$1
trait_type=$2
Rscript run_bmm_m2.R ${SLURM_ARRAY_TASK_ID} ${trait} ${trait_type}
