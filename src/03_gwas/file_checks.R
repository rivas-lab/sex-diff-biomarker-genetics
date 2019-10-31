# file_checks.R
# E Flynn
# 11/17/2017
#
# File checks - checks if a given GWAS file is present.
#  DEPRECATED, will do in bash for sake of speed. 

args <- commandArgs(trailingOnly=TRUE)
trait <- args[1]

GWAS.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/biomarker_gwas/"


fileChecks <- function(file.f, file.m, chr, field){

	missing.list <- data.frame(trait=character(), chr=numeric(), sex=factor(levels=c("zerosex", "onesex")))
	if (!file.exists(file.f)){
        print(sprintf("File missing for c%s trait:%s - female", chr, field))
        missing.list <- rbind(missing.list, list("trait"=field, "chr"=chr, "sex"=1) )
    } else {
	    try.f <- try(read.table(file.f))
	    if (inherits(try.f, "try-error")){
	        print(sprintf("Error loading files for c%s trait:%s - female", chr, field))
	        missing.list <-rbind(missing.list, list("trait"=field, "chr"=chr, "sex"=1) )
	    }
    
    }
    if (!file.exists(file.m)){
        print(sprintf("File missing for c%s trait:%s - male", chr, field))
        missing.list <- rbind(missing.list, list("trait"=field, "chr"=chr, "sex"=2) )
    } else {
		try.m <- try(read.table(file.m))
		if (inherits(try.m, "try-error")){
			print(sprintf("Error loading files for c%s trait:%s - male", chr, field))
			missing.list <- rbind(missing.list, list("trait"=field, "chr"=chr, "sex"=2) )
		}
    }
    #missing.df <- do.call(rbind, missing.list)
	return(missing.list)
}

checkFileQuant <- function(chr, field){
    prefix <- sprintf("%sukb24893_v2.%s", GWAS.FOLDER, field) 
    #file.f <- paste(c(prefix, ".zerosex.PHENO1_c", chr, ".glm.linear.gz"), collapse="")
    file.f <- paste(c(prefix, ".zerosex.", field, "_c", chr, ".glm.linear.gz"), collapse="")
    #file.m <- paste(c(prefix, ".onesex.PHENO1_c", chr, ".glm.linear.gz"), collapse="")
    file.m <- paste(c(prefix, ".onesex.", field, "_c", chr, ".glm.linear.gz"), collapse="")

    my.classes = c("character", "numeric", "character", "character","character", "character",
                   "numeric", "numeric", "numeric", "numeric", "numeric")

    col.labels <- c("CHROM", "POS", "ID", "REF", "ALT1", "TEST", "OBS_CT", 
        "BETA", "SE", "T_STAT", "P")

    checks <- fileChecks(file.f, file.m, chr, field)
    return(checks)
}

checks <- lapply(c(1:22,"X", "XY", "Y"), function(chr) checkFileQuant(chr, trait))
print(do.call(rbind, checks))