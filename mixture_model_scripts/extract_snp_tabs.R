

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
    if ("2" %in% names(cat.count) | "3" %in% names(cat.count)){
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
        filt.f <- list.ds2[[1]]
        filt.m <- list.ds2[[2]]

        snp.tab <- sexSpecSNPtables( filt.f, filt.m, snp.df) 
        f.tab <- annotateSNP(snp.tab$'1')
        m.tab <- annotateSNP(snp.tab$'2')
        write.table(f.tab, file=sprintf("%s/f_spec_%s.txt", DATA.FOLDER2, trait), row.names=FALSE)
        write.table(m.tab, file=sprintf("%s/m_spec_%s.txt", DATA.FOLDER2, trait), row.names=FALSE)
    } else {
        print(sprintf("No sex-specific SNPs for %s", trait))
    }
}

sapply(list.traits, extractSNP)