#!/bin/bash
#SBATCH --job-name=snpnet
#SBATCH --output=logs/snpnet.%A.out
#SBATCH  --error=logs/snpnet.%A.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=150000
#SBATCH --time=1-10:00:00
#SBATCH -p mrivas
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

pop_name=$1
phenotype_name="Testosterone"

snpnet_dir="@@@@@/software/snpnet"
family="gaussian"
geno_dir="@@@@@/geno/array_combined"

out_dir_root="@@@@@/sex-div-analysis/06_snpnet/v1/${pop_name}"
phe_file="@@@@@/sex-div-analysis/06_snpnet/phe_data/${pop_name}.phe"
covariates="None"

if [ ! -d ${out_dir_root} ] ; then mkdir -p ${out_dir_root} ; fi

echo "[$0 $(date +%Y%m%d-%H%M%S)] [start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2
bash @@@@@/helper/snpnet_wrapper.sh ${snpnet_dir} ${phenotype_name} ${family} ${geno_dir} ${out_dir_root} ${phe_file} ${covariates} ${cores} ${mem}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2

