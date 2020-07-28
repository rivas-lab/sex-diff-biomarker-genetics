#!/bin/bash
# run LD for a specific chromosome

N=$1

if [ $# -lt 6 ] ; then memory=64000 ; else memory=$6 ; fi
if [ $# -lt 7 ] ; then threads=8 ; else threads=$7 ; fi

#fill in
target=ukbb_data/
outDir=data/chr_qc/


bed=$target/cal/pgen/ukb24983_cal_chr${N}_v2.bed
bim=$target/snp/snp_download/ukb_snp_chr${N}_v2.bim
fam=$target/fam/ukb2498_cal_v2_s488370_withoutBatch.fam
covar=$target/sqc/ukb24983_GWAS_covar.phe
navars=$target/cal/gwas/na.vars.list # variants only on one array

out=$outDir/ld_out${N}.txt

# PLINK
plink2 --indep-pairwise 100 10 0.1 --memory $memory --threads $threads --bed $bed --bim $bim --fam $fam  --out $out

