# Applying `snpnet` to fit a PRS model

We applied [snpnet](https://github.com/rivas-lab/snpnet) to construct PRS models.

We have improved the snpnet wrapper script and the usage shown in this repository does NOT reflect the most recent syntax/usage of the snpnet wrapper script. For example, we call `snpnet_wrapper.sh` in most of the scripts and that is [now part of the snpnet repository](https://github.com/rivas-lab/snpnet/blob/bd9c27c920ead06b811ce3d69197d1c8967ac1c2/helpers/snpnet_wrapper.sh). Please check the [snpnet](https://github.com/rivas-lab/snpnet) reporitory if you're interested in using this software.

- v1: the original run (the results are NOT used in the final write-up)
  - [`snpnet-v1.20191014.sbatch.sh`](snpnet-v1.20191014.sbatch.sh): a SLURM job script to call `snpnet`
  - [`snpnet-v1.20191014.post-process.sh`](snpnet-v1.20191014.post-process.sh): this script computes PRS for all individuals (including the one in test set) from the BETAs (coefficients) from snpnet.
  - we looped over those scripts for the four population groups: `zerosex`, `onesex`, `post_meno`, and `pre_meno`.
- v2: We updated the list of covariates. We originally passed the covariate values to `snpnet`
  - The job submission script and the post-processing scripts are: [`snpnet-v2.20191026.sbatch.sh`](snpnet-v2.20191026.sbatch.sh) and [`snpnet-v2.20191026.post-process.sh`](snpnet-v2.20191026.post-process.sh)
- v2.1: We decided to pass the residualized phenotypes to `snpnet`.
  - The job submission script and the post-processing scripts are: [`pnet-v2.1.20191026.sbatch.sh`](pnet-v2.1.20191026.sbatch.sh) and [`snpnet-v2.1.20191026.post-process.sh`](snpnet-v2.1.20191026.post-process.sh).
- v2.1-basic-covars: An additional run focusing on a minimal set of covariates (age, [sex,] PC1-PC10).
  - [`snpnet-v2.1-basic-covars.sh`](snpnet-v2.1-basic-covars.sh): this is the job submission script. The post-processing comptuation is now part of the pipeline so there is no need for an extra script.
  - [`snpnet-v2.1-no-covars.sh`](snpnet-v2.1-no-covars.sh): we also experimented with a run without covariates. We did NOT use this results in the final results.
