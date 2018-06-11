#!/bin/bash
#SBATCH --job-name=extract_h
#SBATCH --output=logs/extract_h_%A.out
#SBATCH --error=logs/extract_h_%A.err
#SBATCH --time=0:20:00
#SBATCH -p rbaltman
#SBATCH --mem=4000

Rscript extractUpdatedHerit.R $1