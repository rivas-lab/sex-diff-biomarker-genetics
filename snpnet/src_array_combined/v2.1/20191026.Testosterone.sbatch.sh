#!/bin/bash
#SBATCH --job-name=snpnet
#SBATCH --output=logs/snpnet.%A.out
#SBATCH  --error=logs/snpnet.%A.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=150000
#SBATCH --time=4-0:00:00
#SBATCH -p mrivas
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

pop_name=$1

phenotype_name="Testosterone_residuals"
snpnet_dir="/oak/stanford/groups/mrivas/software/snpnet"
family="gaussian"
#geno_dir="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/split"
geno_dir="/scratch/users/ytanigaw/tmp/snpnet/geno/array_combined"

out_dir_root="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1/${pop_name}"
phe_file="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/phe_data/v2/${pop_name}.phe"
covariates="None"

if [ ! -d ${out_dir_root} ] ; then mkdir -p ${out_dir_root} ; fi

echo "[$0 $(date +%Y%m%d-%H%M%S)] [start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2
bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/PRS/helper/snpnet_wrapper.sh ${snpnet_dir} ${phenotype_name} ${family} ${geno_dir} ${out_dir_root} ${phe_file} ${covariates} ${cores} ${mem}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2

