# map_phe_ids.R
# E Flynn
# 9/11/2017
# 
# Code for mapping phenotype IDs from first to second application.

args = commandArgs(trailingOnly=TRUE)

MAPPING_FILE <- '/scratch/PI/mrivas/ukbb/24983/phe/ukb24983_16698_mapping.tsv' # columns: 24983, 16698

PHE_OUT_DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phefiles'
PHE_IN_DIR <- '/share/PI/mrivas/data/ukbb/phefiles'


# Load the mapping data, convert to a map
map.data <- read.table(MAPPING_FILE, header=FALSE, colClasses="character") # casting to character to avoid ambiguity
colnames(map.data) <- c("curr", "old")
old.to.new <- split(map.data$curr, map.data$old)


mapPheIDs <- function(phe_file, phe_id){
	phe.data <- read.table(sprintf('%s/%s', PHE_IN_DIR, phe_file), header=FALSE, 
		colClasses=c("character", "character", "numeric"))
	colnames(phe.data) <- c("ID", "ID2", "phe")

	print(table(phe.data$ID %in% names(old.to.new)))

	phe.data.filt <- phe.data[phe.data$ID %in% names(old.to.new),]
	phe.data.filt$ID.new <- old.to.new[phe.data.filt$ID] 
	print(summary(phe.data.filt$phe))
	print(table(phe.data.filt$phe <0))
	phe.data.filt <- phe.data.filt[phe.data.filt$phe >0, ] ### remove missing
	phe.data.new <- cbind(phe.data.filt$ID.new, phe.data.filt$ID.new, phe.data.filt$phe)

	write.table(phe.data.new, file=sprintf("%s/%s.phe", PHE_OUT_DIR, phe_id), col.names=FALSE, row.names=FALSE, quote=FALSE)	
}

# for now - using for 
mapPheIDs(args[[1]], args[[2]])







