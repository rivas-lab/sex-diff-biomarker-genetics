#!/bin/bash
#SBATCH --job-name=post_calc
#SBATCH --output=logs/post_calc_%A.out
#SBATCH --error=logs/post_calc_%A.err
#SBATCH --time=0:30:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=8000

ml r-rstan
Rscript extractPosteriorM2.R $1
