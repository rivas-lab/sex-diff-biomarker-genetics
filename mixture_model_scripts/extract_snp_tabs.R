source('model_utils.R')
#source('heritability_utils.R')
source('snp_utils.R')
require('stringr')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"

extractSNPcat <- function(snp.df, df.f, df.m, category, trait){
    comp4 <- snp.df[which(snp.df$category==category),c("p1", "p2", "p3", "p4", "SNP")]
    if (length(comp4$SNP) > 0){
            both.snps <- cbind(df.f[df.f$SNP %in% comp4$SNP ,c("SNP", "CHR", "BP", "BETA","SE", "P")], 
             df.m[df.m$SNP %in% comp4$SNP,c("BETA","SE", "P")])
            colnames(both.snps) <- c("SNP", "CHR", "BP", "B_f", "SE_f", "p_f", "B_m", "SE_m", "p_m")
            both.snp.df <- both.snps[,c("SNP", "CHR", "BP", "B_f", "B_m", "SE_f", "SE_m", "p_m","p_f")] 
            both.snp.df <- merge(both.snp.df, comp4, by="SNP")       
    both.snp.df2 <- annotateSNP(both.snp.df)

    write.table(both.snp.df2, file=sprintf("%s/biomarker/m2/snps%s_%s.txt", DATA.FOLDER, category, trait), row.names=FALSE)

    }   


}

extractSNP <- function(trait){

	 res = tryCatch({
snp.df <- read.delim(sprintf("%s/biomarker/m2/snp_table_%s.txt", DATA.FOLDER, trait), header=TRUE, sep=" ")
print(trait)
cat.count <- table(snp.df$category)
print(cat.count)
    if ("2" %in% names(cat.count) | "3" %in% names(cat.count) | "4" %in% names(cat.count)){
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
    print("Extracting")
    extractSNPcat(snp.df, df.f, df.m, 2, trait)
    extractSNPcat(snp.df, df.f, df.m, 3, trait)
    extractSNPcat(snp.df, df.f, df.m, 4, trait)
    } else {
        print(sprintf("No sex-specific SNPs for %s", trait))
    }
}, error = function(err) {
   print(sprintf("Error for %s", trait))
   print(err)
})

}

args <- commandArgs(trailingOnly=TRUE)
trait.idx <- as.numeric(args[1])
biomarker <- read.table("/scratch/PI/mrivas/users/erflynn/sex_div_gwas/phe_extraction/list_biomarker.txt", header=TRUE, stringsAsFactors=FALSE)

list.traits <- sapply(biomarker$name, function(x) {y <- str_trim(x); z <- gsub(" |-", "_", y); z})
#exists <- sapply(list.traits, function(trait) file.exists(sprintf("%s/summary_dat_%s_2_.txt",DATA.FOLDER,trait)))
#list.traits[!exists]  # which failed - check+re-run?

list.traits2 <- list.traits #[exists]

NUM.BLOCK <- 6
if (trait.idx == length(list.traits2)/NUM.BLOCK){
   selected.traits <- list.traits2[((trait.idx-1)*NUM.BLOCK):length(list.traits2)]
} else{
   selected.traits <- list.traits2[((trait.idx-1)*NUM.BLOCK):((trait.idx)*NUM.BLOCK)]

}


lapply(selected.traits, extractSNP)