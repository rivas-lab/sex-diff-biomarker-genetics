#!/bin/bash


BIOMARKER_PHE_DIR=phefiles/biomarker_m_v2/
for phe_f in `ls $BIOMARKER_PHE_DIR`; do 
  trait="${phe_f%.*}";
  sbatch src/04_bmm/submit_run_bmm.sh ${trait};
  echo $trait; 
done
