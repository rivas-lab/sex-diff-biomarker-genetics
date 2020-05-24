# run_bmm.R
# E Flynn
# 11/17/2017
#
# Updated code for running BMM model for a trait.
#
# Key updates:
#   - sped up model/data loading by removing log_lik calculation because I was not using this
#   - extract heritability results as part of the analysis instead of post-processing
#   - extract results divided by chromosome (intended for X/XY/autosomal sub-analyses)

source('src/04_bmm/model_utils.R')
source('src/04_bmm/heritability_utils.R')
source('src/04_bmm/snp_utils.R')

require('R.utils')
require('data.table')
require('tidyverse')
require('reshape2')
require('parallel')

DATA.FOLDER <- "data/"
GWAS.DIR <- "data/gwas_0522/"
snps.to.keep <- read.table(sprintf("%s/snp_filt_list_wX_v3.txt", DATA.FOLDER), header=FALSE, 
    colClasses="character")

# PARSE ARGUMENTS
args = commandArgs(trailingOnly=TRUE) 

model <- as.numeric(args[1])
trait <- args[2]

# 3rd argument - number of dimensions
if (length(args > 2)){
    ndim <- as.numeric(args[3])
} else {
    ndim <- 2
}

#4th argument - data folder
if (length(args > 3)){
  out_dir <- args[4]   
} else {
  out_dir <- DATA.FOLDER
}

#5th argument - gwas folder
if (length(args > 4)){
  in_dir <- args[5]   
} else {
  in_dir <- GWAS.DIR
}


#


# 6th argument - suffix (used for downsampled, training, validation)
if (length(args) > 5){
    suffix <- args[6]
} else { 
    suffix <- ""
}


if (ndim==3 & model==2){
    print("Warning: M2 is not implemented for 3D yet.")
    exit()
}



se.cutoff <- 0.2

ndim_to_prefixes <- list("2"=c("zerosex", "onesex"), 
    "3"=c("pre_meno", "post_meno", "onesex"),
    "4"=c("under65_f", "over65_f", "under65_m", "over65_m"))

list.prefixes <- ndim_to_prefixes[[as.character(ndim)]]
if (suffix != ""){
   list.prefixes <- sapply(list.prefixes, function(prefix) paste(prefix, suffix, sep="_"))
}

print(sprintf("Running %s for trait %s with %s dim, %s, with prefixes %s", model, trait, ndim, suffix, paste(list.prefixes, collapse=" ")))



getData <- function(trait){
       # for each trait
    list.ds <- lapply(list.prefixes, function(prefix) {
        dat.1 <- fread(sprintf("%s/ukb24983_v2_hg19.%s_%s.genotyped.glm.linear", in_dir, trait, prefix), data.table=FALSE)
        colnames(dat.1) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P")
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
    dat.file <- sprintf("%s/dat_%s.RData", out_dir, trait)

    if (!file.exists(dat.file)){
        list.ds2 <- getData(trait)
    	dat<- extractDataStanMulti(list.ds2)        
    	save(dat, file=sprintf("%s/dat_%s.RData", out_dir, trait))
        } else{
            load(dat.file)
    }

    return(dat)
}




runM1 <- function(trait){
	# run model 1 for a specified trait
    dat <- loadDat(trait)
    dat$dat$K <- 2
    print("LEARNING PARAMS")
    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    fit1 <- stan(file = "src/04_bmm/models/model1_no_loglik.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
    print("SAVING")
    rm(dat)
    print(fit1, pars=c("Sigma", "pi", "Omegacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    save(fit1, file=sprintf("%s/m1/f_%s.RData", out_dir, trait))

}


runM2 <- function(trait){
    # run model 2 for a specified trait
    
    dat <- loadDat(trait)
    dat$dat$K <- 4

    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    fit2 <- stan(file = "tmp_models_tfr/m2_2.stan",  
            data = dat$dat,    
            chains = 4, warmup =1000, iter = 2000, cores = 4, refresh = 200)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma", "lp__"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/m2/f_m2_%s.RData", out_dir, trait))
}

runM2_alt <- function(trait){
    # run model 2 for a specified trait

    dat <- loadDat(trait)
    dat$dat$K <- 2

    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    fit2 <- stan(file = "src/04_bmm/models/model2_alt.stan",  
            data = dat$dat,    
            chains = 4, warmup =200, iter = 600, cores = 4, refresh = 200)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/m2/f_alt_m2_%s.RData", out_dir, trait))
}


extractDataM1 <- function(trait){

    print("EXTRACTING")
    # load the data and fit to extract info about the run
    load(file=sprintf("%s/dat_%s.RData", out_dir, trait))
    load(file=sprintf("%s/m1/f_%s.RData", out_dir, trait))
    
    m1.pi <- getPi(fit1)
    m1.Sigma <- getSigmaMulti(fit1, ndim)
    rg <- getRgMulti(fit1, ndim)
    rg.c <- getRgConfMulti(fit1, ndim)

    rm(fit1)

    # assign each SNP to a component, estimate heritability
    #dat <- labelCategories(dat, m1.Sigma, m1.pi) # label the SNPs with heritability
    #h <- overallHeritability(dat, m1.Sigma, m1.pi)
    h <- c(NA, NA)

    # write out the data
    next.row <- data.frame(t(c(trait, dat$dat$N, unlist(m1.pi), unlist(m1.Sigma), unlist(rg), unlist(rg.c$l), unlist(rg.c$u), unlist(h))))

    # TODO - label these better
    write.table(next.row, file=sprintf("%s/m1/summary_dat_%s_%s_%s.txt", out_dir, trait, ndim, suffix), row.names=FALSE, quote=FALSE)
}

calcErrBarsHerit <- function(trait){
    fit.file=sprintf("%s/m1/f_%s.RData", out_dir, trait)
    if (!file.exists(fit.file)){
            df <- data.frame(t(c("NA", trait, "NA", "NA")))
            colnames(df) <- c("value", "trait", "int", "sex")
            return(df)

        }
        load(file=fit.file)
        load(file=sprintf("%s/dat_%s.RData", out_dir, trait))

        # extract all estimate
        list_of_draws <- rstan::extract(fit1)
        pi.draws <- list_of_draws$pi
        p <- pi.draws
        s.draws <- list_of_draws$Sigma
        Sigma <- s.draws
	
	# extract lower + upper pi
        ordered.p <- p[order(p[,2]),] # ordering p by the non-null component 
        p.lower <- ordered.p[0.025*nrow(ordered.p),]
        p.upper <- ordered.p[0.975*nrow(ordered.p),]
        p.center <- ordered.p[0.50*nrow(ordered.p),]
    
        # extract lower + upper sigma
        ordered.S <- Sigma[order(Sigma[,1,1]),,]
        s.upper <- ordered.S[0.975*dim(Sigma)[1],,]
        s.lower <- ordered.S[0.025*dim(Sigma)[1],,]
        s.center <- ordered.S[0.50*dim(Sigma)[1],,]

        # recalculate SNP membership
        dat2 <- dat
        dat2$categories <- NULL
        dat.u <- labelCategories(dat2, s.upper, p.upper)
        dat.l <- labelCategories(dat2, s.lower, p.lower)
        dat.c <- labelCategories(dat2, s.center, p.center)

        h.up <- overallHeritability(dat.u, s.upper, p.upper)
        h.low <- overallHeritability(dat.l, s.lower, p.lower)
        h.center <- overallHeritability(dat.c, s.center, p.center)

        res <- list("up"=h.up, "low"=h.low, "center"=h.center)

        # reformat into data frame

     my.df <- cbind(t(as.data.frame(res)), trait)
        my.df2 <- data.frame(cbind(my.df, rownames(my.df)))

    colnames(my.df2) <- c("hf", "hm", "trait", "int")
    my.df3 <- melt(my.df2, id.vars=c("trait", "int"), variable.name="sex")
         rownames(my.df3) <- NULL
	write.table(my.df3, file=sprintf("%s/m1/h_err_%s_%s_%s.txt", out_dir, trait, ndim, suffix), row.names=FALSE, quote=FALSE, sep="\t")
    return(my.df3)
    }


### post-processing
extractDataM2 <- function(trait){
    print("Extracting")
    print(trait)

    load(file=sprintf("%s/dat_%s.RData", out_dir, trait)) 
    load(file=sprintf("%s/m2/f_m2_%s.RData", out_dir, trait))

    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    Sigma <- getSigma(fit2)
 
    # assign each SNP to a category
    posterior.df <- posteriorSNPtable(dat, fit2)
    write.table(posterior.df, file=sprintf("%s/m2/snp_table_%s.txt", out_dir, trait), row.names=FALSE, quote=FALSE)
    print("Posterior table generated")

    # remove large files from the workspace
    rm(fit2)

    cat.count <- table(posterior.df$category)
    print(cat.count)

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
    write.table(non.null.snps, file=sprintf("%s/m2/sig_snps%s_.txt", out_dir,trait), row.names=FALSE)
    annot.snp <- annotateSNP(non.null.snps)
     
    sapply(unique(annot.snp$category), function(category) {
        annot.snp.cat <- annot.snp %>% dplyr::filter(category==category)
        write.table(annot.snp.cat, file=sprintf("%s/m2/snps%s_%s.txt", out_dir, category, trait), row.names=FALSE)
         } )
}




if (model==1){
   runM1(trait)
   extractDataM1(trait)
#   calcErrBarsHerit(trait)
}
if (model==2){
    runM2(trait)
   extractDataM2(trait)
}
if (model==3){
   runM2_alt(trait)
}

