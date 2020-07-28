#!/bin/bash
#
# 04_submit_bmm.sh
# E Flynn
#
# Wrapper to submit the mixture model scripts for all traits.

BIOMARKER_PHE_DIR=phefiles/biomarker_m/
for phe_f in `ls $BIOMARKER_PHE_DIR`; do 
  trait="${phe_f%.*}";
  sbatch src/04_bmm/submit_run_bmm.sh ${trait};
  echo $trait; 
done
