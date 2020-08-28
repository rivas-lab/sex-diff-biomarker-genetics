# figshare submission

We used [`export_snpnet_BETA.R`](export_snpnet_BETA.R) to export the PRS coefficients and released it on figshare.

## figshare

- DOI: https://doi.org/10.6084/m9.figshare.12793490

## draft of the data descriptor on figshare

### title

The snpnet polygenic risk score coefficients for Testosterone levels described in 'Sex specific genetic effects across biomarkers'

### description

The dataset contains the coefficients of the polygenic risk scores for Testosterone levels described in the following pre-print:

E. Flynn, Y. Tanigawa, F. Rodriguez, R. B. Altman, N. Sinnott-Armstrong, M. A. Rivas, Sex-specific genetic effects across biomarkers. bioRxiv, 837021 (2019). doi:10.1101/837021

We provide 3 files corresponding to the polygenic risk score models training on the following set of individuals:

- "snpnet.BETAs.Testosterone.combined.tsv.gz": A PRS model trained on both male and female individuals
- "snpnet.BETAs.Testosterone.female-specific.tsv.gz": A PRS model trained on female individuals
- "snpnet.BETAs.Testosterone.male-specific.tsv.gz": A PRS model trained on male individuals

Those files are compressed tab-delimited table files, each of which contains the coefficients (weights) of the polygenic risk score and have the following columns:

- CHROM: the chromosome
- POS: the position
- ID: the variant identifier
- REF: the reference allele
- ALT: the alternate allele
- BETA: the coefficients (weights) of the PRS

Note that we used GRCh37/hg19 genome reference in the analysis and the BETA is always reported for the alternate allele.

We used the BASIL algorithm implemented in R snpnet package, which is described in another preprint:

J. Qian, et al, A Fast and Flexible Algorithm for Solving the Lasso in Large-scale and Ultrahigh-dimensional Problems. bioRxiv, 630079 (2019). doi:10.1101/630079
