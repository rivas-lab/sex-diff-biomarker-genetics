#!/bin/bash
#SBATCH --job-name=bio_v
#SBATCH --output=bio_v.out
#SBATCH --error=bio_v.err
#SBATCH --time=1:00:00 
#SBATCH --partition=rbaltman
#SBATCH --mem=50000


Rscript 02_extract_biomarker_visit.R
