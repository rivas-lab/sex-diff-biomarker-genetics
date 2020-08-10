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
# combined zerosex onesex

if [ $pop_name == "combined" ] ; then
    sex_covars="sex,"
else
    sex_covars=""
fi

phenotype_name="Testosterone"
snpnet_dir="@@@@@/software/snpnet"
family="gaussian"
geno_dir="@@@@@/array_combined"

out_dir_root="@@@@@/sex-div-analysis/06_snpnet/v2/${pop_name}"
phe_file="@@@@@/sex-div-analysis/06_snpnet/phe_data/v2/${pop_name}.phe"
common_covars="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,PC30,PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40,UrineSampleMinutes,DrawTime,DilutionFactorTimeZero,FastingTime,f.53.0.0_2006,f.53.0.0_2007.04,f.53.0.0_2007.05,f.53.0.0_2007.06,f.53.0.0_2007.07,f.53.0.0_2007.08,f.53.0.0_2007.09,f.53.0.0_2007.10,f.53.0.0_2007.11,f.53.0.0_2007.12,f.53.0.0_2008.01,f.53.0.0_2008.02,f.53.0.0_2008.03,f.53.0.0_2008.04,f.53.0.0_2008.05,f.53.0.0_2008.06,f.53.0.0_2008.07,f.53.0.0_2008.08,f.53.0.0_2008.09,f.53.0.0_2008.10,f.53.0.0_2008.11,f.53.0.0_2008.12,f.53.0.0_2009.01,f.53.0.0_2009.02,f.53.0.0_2009.03,f.53.0.0_2009.04,f.53.0.0_2009.05,f.53.0.0_2009.06,f.53.0.0_2009.07,f.53.0.0_2009.08,f.53.0.0_2009.09,f.53.0.0_2009.10,f.53.0.0_2009.11,f.53.0.0_2009.12,f.53.0.0_2010.01,f.53.0.0_2010.02,f.53.0.0_2010.03,f.53.0.0_2010.04,f.53.0.0_2010.05,f.53.0.0_2010.06,f.53.0.0_2010.07,f.53.0.0_2010.0810,f.54.0.0_10003,f.54.0.0_11001,f.54.0.0_11002,f.54.0.0_11003,f.54.0.0_11004,f.54.0.0_11005,f.54.0.0_11006,f.54.0.0_11007,f.54.0.0_11008,f.54.0.0_11009,f.54.0.0_11010,f.54.0.0_11011,f.54.0.0_11012,f.54.0.0_11013,f.54.0.0_11014,f.54.0.0_11016,f.54.0.0_11017,f.54.0.0_11018,f.54.0.0_11020,f.54.0.0_11021,f.54.0.0_11022,f.54.0.0_11023,ageIndicator,ageBin,Batch,Ethnicity,ageBin_FastingTime"

covariates="${sex_covars}${common_covars}"

if [ ! -d ${out_dir_root} ] ; then mkdir -p ${out_dir_root} ; fi

echo "[$0 $(date +%Y%m%d-%H%M%S)] [start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2
bash @@@@@/snpnet_wrapper.sh ${snpnet_dir} ${phenotype_name} ${family} ${geno_dir} ${out_dir_root} ${phe_file} ${covariates} ${cores} ${mem}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; pop_name = ${pop_name}; phenotype = ${phenotype_name}" >&2
