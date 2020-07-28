#!/bin/bash
#SBATCH --job-name=ld
#SBATCH --output=logs/ld_%A_%a.out
#SBATCH --error=logs/ld_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=10000

CHR=$1

./src/02_variant_filtering/02a_run_ld_pairwise.sh $CHR



