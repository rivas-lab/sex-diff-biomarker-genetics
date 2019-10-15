#!/bin/bash
#SBATCH --job-name=run_m1_%A
#SBATCH --output=logs/run_m1_%A_%a.out
#SBATCH --error=logs/run_m1_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=20:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=10000
# usage:
#    sbatch submit_run_bmm.sh <trait> <ndim> <downsampled>

ml r-rstan

trait=$1
ndim=$2

Rscript run_bmm_m1.R ${trait} ${ndim}
