
# Instructions for running mixture model scripts.

### E Flynn
### Updated: 6/10/2019


Once, before running, it is important to run SNP filtering based on Chris' filters and LD structure. 

```Rscript filter_snp.R```

This generates a list of SNPs to keep (~360k) in the data directory. 

Next, we can run m1, using the script:

```submit_run_bmm.sh```

For now, edit the file directly to contain the appropriate arguments for the trait of interest. This will be updated.
This submit script is a wrapper for run_bmm_m1.R. 


To run model 2 (the four component model):

```sbatch run_bmm_m2.sh $trait "quant" ```

This runs the script run_bmm_m2.R. This script has three main functions:
loadDat, runM2, and extractData. The script also uses functions from model_utils.R and SNP_utils.R.

The model used is model2.stan. This model does -NOT- use a variance covariance matrix for the fourth component (model2_v2.stan does, but there were some identifiability issues). This may be a cause of some of the problems.