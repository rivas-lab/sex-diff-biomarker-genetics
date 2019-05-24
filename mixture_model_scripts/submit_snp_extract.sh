#!/bin/bash
#SBATCH --job-name=h_snp
#SBATCH --output=logs/h_snp_%A_%a.out
#SBATCH --error=logs/h_snp_%A_%a.err
#SBATCH --time=1:30:00
#SBATCH -p rbaltman
#SBATCH --mem=10000

ml r-rstan
Rscript extract_snp_tabs.R

