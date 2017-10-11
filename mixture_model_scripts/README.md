
# Instructions for running mixture model scripts.

### E Flynn
### Updated: 10/10/2017


Once, before running, it is important to run SNP filtering based on Chris' filters and LD structure. 

```Rscript filter_snp.R```

This generates a list of SNPs to keep (~360k) in the data directory. 

Next, we can run the models, using the script:

```submit_run_bmm.sh```

For now, edit the file directly to contain the appropriate arguments for the trait of interest. This will be updated.
This submit script is a wrapper for run_bmm.R. 