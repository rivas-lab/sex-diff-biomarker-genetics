#!/bin/bash
set -beEuo pipefail


snpnet_post () {
    local split=$1
    local phe=$2

    local covars="None"
    local data_dir_root="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out"
    local helper_dir="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/PRS/helper"
    local geno="/oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv"

    bash ${helper_dir}/export_betas.sh \
        ${data_dir_root}/${split} ${phe} ${covars}

    bash ${helper_dir}/plink_score.sh \
        ${data_dir_root}/${split} ${phe} ${geno}
}

split=$1
phe="Testosterone"

snpnet_post $split $phe


