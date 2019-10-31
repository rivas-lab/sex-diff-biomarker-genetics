#!/bin/bash
#SBATCH --job-name=post_calc
#SBATCH --output=logs/post_calc.out
#SBATCH --error=logs/post_calc.err
#SBATCH --time=1:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=10000

ml r-rstan
Rscript extractPosteriorM1.R
