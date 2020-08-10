# get_xy_chr_stats.sh
# E Flynn
# 11/17/2017
#
# Code used for getting X, Y, XY, MT chromosome stats. 
# Gathers LD, frequency, and missingness info using plink. 


N=${1} # the chromosome (X, Y, XY, MT)



# updated to match other chromosomes
target=ukbb/24983

bed=$target/cal/pgen/ukb24983_cal_chr${N}_v2.bed
bim=$target/snp/snp_download/ukb_snp_chr${N}_v2.bim
fam=$target/fam/ukb2498_cal_v2_s488370_withoutBatch.fam
out=/scratch/users/erflynn/bmm_project/data/chr_qc/chr${N}

#plink --bed ${bed} --bim ${bim} --fam ${fam} --indep 50 5 2 --out ${out}_ld #submitted via sumbit_ld.sh
plink2 --bed ${bed} --bim ${bim} --fam ${fam} --hardy --out ${out}_hardy
plink2 --bed ${bed} --bim ${bim} --fam ${fam} --freq --out ${out}.f
plink2 --bed ${bed} --bim ${bim} --fam ${fam} --missing --out ${out}.m

