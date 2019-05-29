

require('rstan')
require('ggplot2')
require('tidyverse')
require('reshape2')
source('model_utils.R')
source('snp_utils.R')

DATA.FOLDER2 <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/biomarker/m2"

biomarker <- read.table("/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phe_extraction/list_biomarker.txt", header=TRUE, stringsAsFactors=FALSE)

list.traits <- biomarker$Field
#m2exists <- sapply(list.traits, function(trait) file.exists(sprintf("%s/snp_table_%s.txt",DATA.FOLDER2,trait)))
#list.traits2 <- list.traits[!m2exists]
print(list.traits)

extractSNP <- function(trait){
     snp.df <- read.table(sprintf("%s/snp_table_%s.txt",DATA.FOLDER2, trait), header=TRUE, stringsAsFactors=FALSE)
    cat.count <- table(snp.df$category)
    
        print(trait)
        print(cat.count)
        list.prefixes <- c("zerosex", "onesex")
        chrs <- c(1:22)
        list.ds <- lapply(list.prefixes, function(prefix) {
		all.dat <- do.call(rbind, lapply(chrs, function(chr) { getFile(prefix, chr, trait)}));
                colnames(all.dat)[1:3] <- c("CHR", "BP", "SNP");
				       return(all.dat)
            })

        list.ds2 <- extractOverlappingRows(list.ds)
	df.f <- list.ds2[[1]]
	df.m <- list.ds2[[2]]

    comp4 <- snp.df[which(snp.df$category==4),c("p1", "p2", "p3", "p4", "SNP")]
    comp4snps <- comp4$SNP
    if (length(comp4snps) > 0){
            both.snps <- cbind(df.f[df.f$SNP %in% comp4snps ,c("SNP", "CHR", "BP", "BETA","SE", "P")], 
             df.m[df.m$SNP %in% comp4snps,c("BETA","SE", "P")])
            colnames(both.snps) <- c("SNP", "CHR", "BP", "B_f", "SE_f", "p_f", "B_m", "SE_m", "p_m")
            both.snp.df <- both.snps[,c("SNP", "CHR", "BP", "B_f", "B_m", "SE_f", "SE_m", "p_m","p_f")] 
            both.snp.df <- merge(both.snp.df, comp4, by="SNP")       
    }   
    both.snp.df2 <- annotateSNP(both.snp.df)

    write.table(both.snp.df2, file=sprintf("%s/both_%s.txt", DATA.FOLDER2, trait), row.names=FALSE)
}

sapply(list.traits, extractSNP)