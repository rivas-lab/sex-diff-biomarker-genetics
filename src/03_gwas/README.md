# Code for sex-divided GWAS runs
### E Flynn
### 9/11/2017



Run this once to gather sex labels and filter:

```Rscript retrieve_sex_labels.R```


Here are the scripts I ran. For now, I am just testing using BMI.





Then, to run GWAS:

```sbatch --array=1-25 submit_gwas_run.sh zerosex```

```sbatch --array=1-25 submit_gwas_run.sh onesex```

```This runs the script gwas_qt_sex_div.sh```

UPDATED
```sbatch --array=1-22 submit_gwas_run.sh <sex> <phe>```


Binary phenotypes

```Rscript map_phe_ids.R rohit/RH107.phe RH107```

```Rscript map_phe_ids.R rohit/RH160.phe RH160```


Then, to run GWAS, do same as above, but use:

```submit_gwas_run_bin.sh```


[deprecated: `map_phe_ids.R`: for mapping btw different versions]