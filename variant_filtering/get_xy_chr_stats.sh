# get_xy_chr_stats.sh
# E Flynn
# 11/17/2017
#
# Code used for getting X, Y, XY, MT chromosome stats. 
# Gathers LD, frequency, and missingness info using plink. 


N=${1} # the chromosome (X, Y, XY, MT)


export PATH=/oak/stanford/groups/mrivas/software/plink2/20170904/:$PATH

#export PATH=/oak/stanford/groups/mrivas/software/plink1.9/plink1.90b4.6/:$PATH # using plink1.9 for indep
# updated to match other chromosomes
target=/oak/stanford/groups/mrivas/ukbb/24983

bed=$target/cal/pgen/ukb24983_cal_chr${N}_v2.bed
bim=$target/snp/snp_download/ukb_snp_chr${N}_v2.bim
fam=$target/fam/ukb2498_cal_v2_s488370_withoutBatch.fam
out=/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/chr_qc/chr${N}

#plink2 --bed ${bed} --bim ${bim} --fam ${fam} --indep-pairwise 100 10 0.1 --out ${out}_ld #submitted via sumbit_ld.sh
plink2 --bed ${bed} --bim ${bim} --fam ${fam} --freq --out ${out}.f
plink2 --bed ${bed} --bim ${bim} --fam ${fam} --missing --out ${out}.m
