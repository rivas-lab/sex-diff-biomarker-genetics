#!/bin/bash
set -beEuo pipefail

pvar_f="@@@@@@@@@@@@@@@@@@@@/ukb24983_cal_hla_cnv.pvar.zst"
data_d="@@@@@@@@@@@@@@@@@@@@"

ml load R/3.6 gcc htslib

sex_stratification_strs=(
combined
male-specific
female-specific
)

for sex_stratification_str in "${sex_stratification_strs[@]}" ; do

in_f="${data_d}/${sex_stratification_str}/Testosterone_residuals/results/betas.tsv"
out_f="${data_d}/snpnet.BETAs.Testosterone.${sex_stratification_str}.tsv"

Rscript /dev/stdin ${in_f} ${out_f} ${pvar_f} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
args <- commandArgs(trailingOnly=TRUE)

in_f   <- args[1]
out_f  <- args[2]
pvar_f <- args[3]

fread(cmd=paste('zstdcat', pvar_f)) %>%
rename('CHROM'='#CHROM') %>%
right_join(
    fread(in_f),
    by=c('ID'='ID', 'ALT'='A1')
) %>%
arrange(CHROM, POS) %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)
EOF

bgzip -l9 -f ${out_f}

echo ${out_f}.gz

done
