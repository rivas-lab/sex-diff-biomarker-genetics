#!/bin/bash
#SBATCH --job-name=run_m2
#SBATCH --output=logs/run_m2_%a_%A.out
#SBATCH --error=logs/run_m2_%a_%A.err
#SBATCH --time=20:00:00
#SBATCH --nodes=4
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=10000
# usage:
#    sbatch run_bmm_m2.sh <trait> <ndim> 


ml r-rstan

trait=$1
ndim=$2

Rscript run_bmm_m2.R ${trait} ${ndim}



