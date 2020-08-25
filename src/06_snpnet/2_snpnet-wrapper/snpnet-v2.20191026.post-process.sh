#!/bin/bash
set -beEuo pipefail

snpnet_post () {
    local split=$1
    local phe=$2

    local covars="None"
    local data_dir_root="@@@@@/sex-div-analysis/06_snpnet/v2"
    local helper_dir="@@@@@/helper"
    local geno="@@@@@/pgen/ukb24983_cal_hla_cnv"


    bash ${helper_dir}/export_betas.sh \
        ${data_dir_root}/${split} ${phe} ${covars}

    bash ${helper_dir}/plink_score.sh \
        ${data_dir_root}/${split} ${phe} ${geno}
}

split=$1
phe="Testosterone"

snpnet_post $split $phe
