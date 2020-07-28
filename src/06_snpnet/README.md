# Polygenic risk prediction of Testosterone level

- `phe_data`: scripts & notebooks to prepare input file (phenotype file) for `snpnet`
- `src_array_combined`: A script to call `snpnet`
- `out`: Outputs from `snpnet`
- `plots`: Notebook to generate plots

In each directory, we have multiple versions
- `v1`: An old version prepared during ASHG2019. This version used a different set of covariates when comparing the sex-specific model and the combined model.
- `v2`: A revised version with a fix for the covariates inconsistency issue. We tried to ask `snpnet` to deal with the covariates.
- `v2.1`: the final results, where we passed the residuals from the covariate correction.

## File location

### PRS score files

```
/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1/combined/Testosterone_residuals/results/score/Testosterone_residuals.sscore
/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1/onesex/Testosterone_residuals/results/score/Testosterone_residuals.sscore
/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1/zerosex/Testosterone_residuals/results/score/Testosterone_residuals.sscore
```

### BETAs from snpnet

```
/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1/combined/Testosterone_residuals/results/betas.tsv
/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1/onesex/Testosterone_residuals/results/betas.tsv
/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/sex-div-analysis/snpnet/out/v2.1/zerosex/Testosterone_residuals/results/betas.tsv
```
