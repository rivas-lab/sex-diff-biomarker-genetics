#!/bin/bash
# run_sex_div_gwas_mult_traits_qt.sh
# E Flynn
# 11/17/2017

# Runs sex divided gwas for multiple traits from an input file.
# The input file is formatted as:
#	PHE_ID   UKBB_TRAIT_ID	 PHE_LOC	DESCRIPTION


BASE_DIR=/scratch/PI/mrivas/users/erflynn/sex_div_gwas
PHE_DIR=${BASE_DIR}/phefiles/processed

log_file=${BASE_DIR}/gwas_scripts/logs/gwas_log.txt

echo "Starting GWAS run" > ${log_file}

cd ${BASE_DIR}/gwas_scripts/

## filter phe files
while IFS=$'\t' read trait trait_id path descript; 
do 
	echo $trait, $path >> ${log_file}; 

	## FILTER TO REMOVE NAs
	path_new=${PHE_DIR}/${trait}.phe 
	Rscript filterPhe.R $path $path_new $trait >> ${log_file} 
	## CHECK SEX-BALANCE FOR TRAIT
	Rscript checkSexBalance.R $path_new $trait >> ${log_file}

	sbatch --array=1-26 submit_gwas_run.sh zerosex $trait $path_new ;
	sbatch --array=1-26 submit_gwas_run.sh onesex $trait $path_new ;

done < ${BASE_DIR}/data/multi_run_gwas_input.txt 

#### second set
log_file=${BASE_DIR}/gwas_scripts/logs/gwas_log_set2.txt

echo "Starting GWAS run" > ${log_file}

while IFS=$'\t' read table trait descript; 
do 
	echo $trait >> ${log_file}; 
	path=${BASE_DIR}/phefiles/${trait}.phe

	## FILTER TO REMOVE NAs
	path_new=${PHE_DIR}/${trait}.phe 
	Rscript filterPhe.R $path $path_new $trait >> ${log_file} 
	## CHECK SEX-BALANCE FOR TRAIT
	Rscript checkSexBalance.R $path_new $trait >> ${log_file}

	sbatch --array=1-26 submit_gwas_run.sh zerosex $trait $path_new ;
	sbatch --array=1-26 submit_gwas_run.sh onesex $trait $path_new ;
done < ${BASE_DIR}/data/set_traits2.txt

###########


#### file checks - make sure files are present and not empty
while IFS=$'\t' read trait trait_id path descript; 
do
	echo $trait;
	./file_checks.sh ${trait} >> ${log_file}
	cp ${BASE_DIR}/results_test/*${trait}*.gz ${BASE_DIR}/results/

done < ${BASE_DIR}/data/multi_run_gwas_input.txt

echo "" > file_checks.out

while IFS=$'\t' read table trait descript; 
do
	echo $trait >> file_checks.out
	./file_checks.sh ${trait} >> file_checks.out
	cp ${BASE_DIR}/results_test/*${trait}*.gz ${BASE_DIR}/results/

done < ${BASE_DIR}/data/set_traits2.txt

###### run BMM


while IFS=$'\t' read trait trait_id path descript; 
do
	echo $trait;
	sbatch submit_run_bmm.sh ${trait};
done < ${BASE_DIR}/data/multi_run_gwas_input.txt

#### TODO - clean up files

# set 2
while IFS=$'\t' read table trait descript;
do
	echo $trait;
	sbatch submit_run_bmm.sh ${trait};
done < ${BASE_DIR}/data/set_traits2.txt