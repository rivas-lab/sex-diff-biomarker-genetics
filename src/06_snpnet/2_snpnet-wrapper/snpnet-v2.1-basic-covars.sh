#!/bin/bash
#SBATCH --job-name=snpnet
#SBATCH --output=logs/snpnet.%A.out
#SBATCH  --error=logs/snpnet.%A.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=80000
#SBATCH --time=4-0:00:00
#SBATCH -p mrivas
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

pop_name=$1
# combined 
# onesex
# zerosex
if [ $pop_name == "combined" ] ; then
    covariates="sex,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
else
    covariates="age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
fi
phenotype_name="Testosterone"
phe_file="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/phe_data/v2.1-no-covars/Testosterone.phe"
genotype_pfile="/scratch/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv"
results_dir="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1-no-covars/${pop_name}_basic-covars"

if [ ! -d ${results_dir} ] ; then mkdir -p ${results_dir} ; fi

echo "[$0 $(date +%Y%m%d-%H%M%S)] [start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2

# ml load snpnet_yt/0.3.5
ml load snpnet_yt/0.3.6

bash ${snpnet_wrapper} --nCores ${cores} --memory ${mem} \
    --covariates ${covariates} --split_col "split_${pop_name}" \
    --verbose --glmnetPlus --save_computeProduct \
    ${genotype_pfile} ${phe_file} ${phenotype_name} \
    "gaussian" ${results_dir}

echo "[$0 $(date +%Y%m%d-%H%M%S)] [end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2
