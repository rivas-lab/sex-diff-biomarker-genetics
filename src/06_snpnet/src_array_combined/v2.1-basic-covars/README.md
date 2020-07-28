# additional `snpnet` runs with minimum sets of covariates

## results files

`/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1-basic-covars/`

## notes

Reviewer's comment:

```{text}
The part about the sex-specific multivariate polygenic prediction is a bit unclear to me. What was actually used in the model? All 1 million variants, or only the ones identified using SEMM? What are "covariate adjusted phenotypes"? Does this refer to the adjustment using 127 variables? If so, I think it would make a lot more sense to make the prediction model based on the unadjusted testosteron values. In a real world setting, it is highly unlikely that those 127 variables will be available for adjustment. I would also be quite curious if the authors report the RMSE or mean error on the untransformed scale.
```

For this, Yosuke decided to add two sets of `snpnet` runs.

1. No-covariates for the log-transformed Testosterone value
2. Minimum set of covariates (age, PC1-PC10) for the log-transformed Testosterone value.

## scripts

- [`phe_file_prep.ipynb`](phe_file_prep.ipynb): phe file prep.
  - Since our initial application of the snpnet, we changed the input file format and the options for the snpnet wrapper. This notebook prepares phe file in new format. We use 3 different split columns to specify training/validation/test sets for male- and female-specific models as well as the combined model.
- [`snpnet-basic-covars.sh`](snpnet-basic-covars.sh): snpnet SLURM job script for the minimum s et of covariates model. We use age, PC1-PC10 for male- and female-specific models, and age, sex, PC1-PC10 for the combined model.
- [`snpnet-no-covars.sh`](snpnet-no-covars.sh): snpnet SLURM job script for no-covariate model. We use no covariate for the sex-specific models and sex in the combined model.

