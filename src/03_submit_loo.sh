#!/bin/bash

for param_id in $(seq 1 14); do
  sbatch --array=1-10 src/05_loo/run_ll.sh 1 $param_id;
  echo $param_id;
done

for param_id in $(seq 1 8); do
  sbatch --array=1-10 src/05_loo/run_ll.sh 2 $param_id;
  echo $param_id;
done
