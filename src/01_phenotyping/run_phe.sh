#!/bin/bash

ml python/2.7.13
ml py-numpy/1.14.3_py27
ml R/3.6.1

# extract
python2 01_extract_med_fields.py 
python2 01_extract_ss_fields.py
python2 01_extract_age_fields.py

# tidy
Rscript 02_extract_med_dat.R
Rscript 02_extract_ss_dat.R
sbatch 02_submit_extract_biomarker.sh 
