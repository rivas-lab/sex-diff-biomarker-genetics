# file_checks.sh
# E Flynn
# 11/17/2017
#
# File checks - checks if a given GWAS file is present.


TRAIT=$1

GWAS_FOLDER=/scratch/PI/mrivas/users/erflynn/sex_div_gwas/biomarker_gwas/

for CHR in `seq 1 26`
do
	f_path=${GWAS_FOLDER}/ukb24893_v2.${TRAIT}.zerosex.PHENO1_c${CHR}.glm.linear.gz;
	m_path=${GWAS_FOLDER}/ukb24893_v2.${TRAIT}.onesex.PHENO1_c${CHR}.glm.linear.gz;

	if [[ ! -s ${f_path} ]]; then echo "File for F c${CHR} is missing for ${TRAIT}"; fi
	if [[ ! -s ${m_path} ]]; then echo "File for M c${CHR} is missing for ${TRAIT}"; fi
	
done

