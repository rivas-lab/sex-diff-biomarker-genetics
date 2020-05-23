#!/bin/bash
#SBATCH --job-name=run_bmm_%A
#SBATCH --output=logs/run_bmm_%A_%a.out
#SBATCH --error=logs/run_bmm_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=24:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=30000
# usage:
#    sbatch submit_run_bmm.sh <model> <trait> <ndim>  <outdir> <indir>



ml purge
ml load R/3.5.1
ml load gcc/8.1.0

model=$SLURM_ARRAY_TASK_ID
trait=Testosterone
ndim=2
outdir=data/
indir=data/gwas
#suffix=$5
Rscript src/04_bmm/run_bmm_vary_params.R ${model} ${trait} ${ndim} ${outdir} ${indir} #${suffix}
