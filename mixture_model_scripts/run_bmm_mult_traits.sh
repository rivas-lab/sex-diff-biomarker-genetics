
# INPUT = file formatted
# PHE_ID   TRAIT_NAME	 PHE_LOC

# step 0 
# -- read in input
while IFS=$'\t' read trait trait_id path descript; 
do 
echo $trait; 
 sbatch submit_run_bmm.sh $trait;
done < ../data/multi_run_gwas_input.txt


# step 1
# -- run mixture model - 1000 iterations?

# step 2
# -- parse output params, rg, heritability
# -- look at x chromosome contributions
