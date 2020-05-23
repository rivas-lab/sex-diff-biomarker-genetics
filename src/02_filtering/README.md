Code for variant filtering steps.

[Note, as written, all code should be run from the base directory, not within this directory.]

First, grab the variant table: \\
  `Rscript 00_filter_variants.R`

Then get the xy chromosome stats and run LD filtering: \\
  `./01_get_xy_chr_stats.sh ${CHR}` \\
  `sbatch 02_submit_ld.sh ${CHR}`


After this, we put things together in a table: \\
 `Rscript 03_qc_filtering.R` \\
and then filter: \\
 `Rscript 04_qc_filtering_part2.R`