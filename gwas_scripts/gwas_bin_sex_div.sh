#!/bin/bash
#
# Original script: Yosuke
# Updated by E Flynn to analyze sex-divided data
# 9/11/2017
#   Changes:
#    - using input param 2 for sex 
#    - using a "sex_ids.keep" file of IDs to keep, this has 151k individuals removed and contains only IDs for samples from one sex 
#    - removed sex as a covariate from list of covariate names
#    - updated file labeling to include sex label
#
set -beEu -o pipefail

############################## START parse args ########################
if [ $# -lt 5 ] ; then echo "usage: $0 <UKBB application ID> <sample_sex> <phen code> <pheno file> <chr = 1-22, X(23), Y(24), XY(25), MT(26)> [memory = 64000] [threads = 8]" >&2 ; exit 1 ; fi
if [ $# -lt 6 ] ; then memory=64000 ; else memory=$6 ; fi
if [ $# -lt 7 ] ; then threads=8 ; else threads=$7 ; fi

UKB_app_id=$1

target=/scratch/PI/mrivas/ukbb/24983
outDir=/scratch/PI/mrivas/users/erflynn/sex_div_gwas/results ### using this for now, update
if [ ! -d $outDir ] ; then mkdir -p $outDir ; fi

# which sex to analyze
sex=$2  # onesex or zerosex

# phenofile
phen=$3
phenofile=$4

# N: chromosome
case "$5" in
	23) N="X" ;;
	24) N="Y" ;;
	25) N="XY" ;;
	26) N="MT" ;;
	*)  N=$5 ;;
esac

############################## END parse args ##########################
which plink2 >&2
which bgzip >&2
which python >&2
which tabix >&2

flipfix_script=/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb24983wiki/scripts/flipfix_A1A2.py  #$(readlink -e $(dirname $(readlink -e $0))/flipfix_A1A2.py)


# filename
bed=$target/cal/ukb_cal_chr${N}_v2.bed
bim=$target/snp/ukb_snp_chr${N}_v2.bim
fam=$target/fam_24983/ukb2498_cal_v2_s488370_withoutBatch.fam
covar=$target/phe/ukb24983_GWAS_covar.phe
#sampleqc=$target/phe/ukb24983_remove.phe
#reducted=$target/phe/w2498_20170726.phe
sex_ids=/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles/${sex}.keep ### sampleqc, reducted already removed from these

for input_f in $phenofile $bed $bim $fam $covar $sex_ids ; do if [ ! -f $input_f ] ; then
	echo "[$0 $(date +%Y%m%d-%H%M%S)] $input_f is missing!" >&2
	exit 1
fi ; done

out="$outDir/ukb${UKB_app_id}_v2_c${N}.${phen}.${sex}"
out_suffix="PHENO1.glm.logistic.hybrid.gz.tbi"
out_renamed="$outDir/ukb${UKB_app_id}_v2.${phen}.${sex}.PHENO1_c${N}.glm.logistic.hybrid.gz.tbi"

# check if the out file exist
for outfile in ${out}.${out_suffix%.gz.tbi} ${out_renamed} ; do
	if [ -f $outfile ] ; then
		echo "[$0 $(date +%Y%m%d-%H%M%S)] the results file $outfile already exists! skipping" >&2
		exit 0
	fi
done

# PLINK
if [ ! -f ${out}.${out_suffix%.gz.tbi} ] ; then
	echo "[$0 $(date +%Y%m%d-%H%M%S)] Running GWAS with plink2" >&2
	plink2  --memory $memory --threads $threads \
		--bed $bed --bim $bim --fam $fam \
		--keep $sex_ids \
		--out $out \
		--pheno $phenofile \
		--covar $covar --covar-name age Array PC1 PC2 PC3 PC4 \
		--glm firth-fallback
fi

# flip fix and rename
echo "[$0 $(date +%Y%m%d-%H%M%S)] flip fix script" >&2
python $flipfix_script -i ${out}.${out_suffix%.gz.tbi} -c $N -b | awk 'NR==1 || $7 == "ADD"' > ${out_renamed%.gz.tbi}

# bgzip & tabix
if [ -f ${out_renamed%.tbi} ] && [ ! -f ${out_renamed%.gz.tbi} ] ; then
	echo "[$0 $(date +%Y%m%d-%H%M%S)] bgzip was incomplete. deleting ${out_renamed%.tbi} and re-do" >&2
	rm ${out_renamed%.tbi}
fi
echo "[$0 $(date +%Y%m%d-%H%M%S)] bgzip and tabix" >&2
bgzip "${out_renamed%.gz.tbi}"
tabix -b2 -e2 "${out_renamed%.tbi}"