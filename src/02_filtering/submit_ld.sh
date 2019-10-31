#!/bin/bash
#SBATCH --job-name=ld
#SBATCH --output=logs/ld_%A_%a.out
#SBATCH --error=logs/ld_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH -p rbaltman
#SBATCH --mem=10000

CHR=$1

export PATH=/oak/stanford/groups/mrivas/software/plink2/20170904/:$PATH
#export PATH=/home/erflynn/applications/plink2/20180223/:$PATH # location of latest plink version
export PATH=$PATH:$HOME/applications/htslib-1.5/installation/bin

./run_ld_pairwise.sh $CHR



