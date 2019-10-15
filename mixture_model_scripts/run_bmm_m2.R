
source('model_utils.R')
source('snp_utils.R')
require('tidyverse')
require('data.table')
require('parallel')

DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/1015/m2/"
GWAS.DIR <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/test_gwas/"

args <- commandArgs(trailingOnly=TRUE)
trait <- args[1]

# 2nd argument - number of dimensions
if (length(args > 1)){
    ndim <- as.numeric(args[2])
} else {
    ndim <- 2
}



maf.cutoff <- 0.01 
se.cutoff <- 0.2


calcLoglik <- FALSE


ndim_to_prefixes <- list("2"=c("zerosex", "onesex"), 
     "3"=c("pre_meno", "post_meno", "onesex"),
     "4"=c("under65_f", "over65_f", "under65_m", "over65_m"))

list.prefixes <- ndim_to_prefixes[[as.character(ndim)]]


getData <- function(trait){
       # for each trait
    list.ds <- lapply(list.prefixes, function(prefix) {
        dat.1 <- fread(sprintf("%s/ukb24983_v2_hg19.%s_%s.genotyped.glm.linear", GWAS.DIR, trait, prefix), data.table=FALSE)
        
        colnames(dat.1)[1:3] <- c("CHR", "BP", "SNP");
        
        rownames(dat.1) <- dat.1$SNP

        # remove NAs
        dat.2 <- dat.1[!is.na(dat.1$SE),]

        # SE filter
        dat.3 <- dat.2[dat.2$SE < QUANT.SE.CUTOFF,]

        # MAF filter
        dat.4 <- dat.3[dat.3$SNP %in% snps.to.keep$V1,]
        return(dat.4)

    })

    list.ds2 <- extractOverlappingRows(list.ds)

}

loadDat <- function(trait){

    list.ds2 <- getData(trait)
    stan.obj <- extractDataStanMulti(list.ds2)
    return(stan.obj)
}


runM2 <- function(trait){
	# run model 2 for a specified trait

    dat <- loadDat(trait)
    dat$dat$K <- 4
    save(dat, file=sprintf("%s/dat_%s.RData", DATA.FOLDER, trait))
    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    if (calcLoglik==TRUE){
        model.file <- "models/model2_loglik.stan" # this is v2...
    } else {
        model.file <- "models/model2.stan" # v2??
    }
    fit2 <- stan(file = model.file,  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 20)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/f_m2_%s.RData", DATA.FOLDER, trait))
}





### post-processing
extractData <- function(trait){
	print("Extracting")
    print(trait)

	load(file=sprintf("%s/dat_%s.RData", DATA.FOLDER, trait)) 
    load(file=sprintf("%s/f_m2_%s.RData", DATA.FOLDER, trait))

    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    Sigma <- getSigma(fit2)
 
    # assign each SNP to a category
    posterior.df <- posteriorSNPtable(dat, fit2)
    write.table(posterior.df, file=sprintf("%s/snp_table_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
    print("Posterior table generated")

    # remove large files from the workspace
	rm(fit2)

    cat.count <- table(posterior.df$category)

    if (ndim!=2){
        print("warning - extraction has yet to be implemented for larger dimensions.")
    }

    # create a data frame with the pvalues
    p.df <- dat$p 
    se.df <- dat$dat$SE
    beta.df <- dat$dat$B
    snp.df <- data.frame(dat$snp)
    chr.df <- data.frame(dat$chr)
    full.df <- do.call(cbind, list(snp.df, chr.df, beta.df, se.df, p.df))
    colnames(full.df) <- c("SNP", "CHR", "B.f", "B.m", "SE.f", "SE.m", "P.f", "P.m")
    comb.df <- cbind(full.df, posterior.df %>% dplyr::select(p1, p2, p3, p4, category))
    non.null.snps <- comb.df %>% dplyr::filter(category %in% c(2,3,4))

    annot.snp <- annotateSNP(non.null.snps)
     
    sapply(unique(annot.snp$category), function(category) {
        annot.snp.cat <- annot.snp %>% dplyr::filter(category==category)
        write.table(annot.snp.cat, file=sprintf("%s/snps%s_%s.txt", DATA.FOLDER, category, trait), row.names=FALSE)
         } )
}


runM2(trait)
extractData(trait)

