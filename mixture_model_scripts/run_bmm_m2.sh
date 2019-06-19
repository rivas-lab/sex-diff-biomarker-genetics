#!/bin/bash
#SBATCH --job-name=run_m2
#SBATCH --output=logs/run_m2_%a_%A.out
#SBATCH --error=logs/run_m2_%a_%A.err
#SBATCH --time=48:00:00
#SBATCH --partition=mrivas
#SBATCH --mem=50000

ml r-rstan

trait=$1
trait_type=$2

if [ "$#" -ne 3 ]; then
    echo "chr number not specified"
    Rscript run_bmm_m2.R ${SLURM_ARRAY_TASK_ID} ${trait} ${trait_type} ;
    exit
fi

chr=$3
echo "running chr " ${chr}
Rscript run_bmm_m2.R ${SLURM_ARRAY_TASK_ID} ${trait} ${trait_type} ${chr} ;



