#!/bin/bash
#SBATCH --job-name=vary_pars
#SBATCH --output=vary_pars_%A.out
#SBATCH --error=vary_pars_%A.err
#SBATCH --time=48:00:00 
#SBATCH --partition=normal
#SBATCH --mem=10000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erflynn@stanford.edu


trait=$1
Rscript test_m1_vary_pars.R ${trait}
