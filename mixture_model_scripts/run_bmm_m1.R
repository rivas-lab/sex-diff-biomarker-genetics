# run_bmm_m1.R
# E Flynn
# 11/17/2017
#
# Updated code for running BMM model #1 a trait.
#
# Key updates:
#   - sped up model/data loading by removing log_lik calculation because I was not using this
#   - extract heritability results as part of the analysis instead of post-processing
#   - extract results divided by chromosome (intended for X/XY/autosomal sub-analyses)

source('model_utils.R')
source('heritability_utils.R')

require('R.utils')
require('data.table')
require('reshape2')
require('parallel')

DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/1015/m1/"
GWAS.DIR <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/test_gwas/"

# PARSE ARGUMENTS
args = commandArgs(trailingOnly=TRUE) # 1, 2, or 3 

trait <- args[1]

# 2nd argument - number of dimensions
if (length(args > 1)){
    ndim <- as.numeric(args[2])
} else {
    ndim <- 2
}

# 3rd argument - downsampled?
if (length(args) > 2){
    downsampled <- args[3]
} else { 
downsampled <- FALSE 
}


maf.cutoff <- 0.01
se.cutoff <- 0.2

ndim_to_prefixes <- list("2"=c("zerosex", "onesex"), 
    "3"=c("pre_meno", "post_meno", "onesex"),
    "4"=c("under65_f", "over65_f", "under65_m", "over65_m"))

list.prefixes <- ndim_to_prefixes[[as.character(ndim)]]

if (downsampled==TRUE){ 
    list.prefixes <- sapply(list.prefixes, function(prefix) paste(prefix, "d", sep="_"))
}

downsampled_str <- ifelse(downsampled, "downsampled", "")
print(sprintf("Running M1 for trait %s with %s dim, %s, with prefixes %s", trait, ndim, downsampled_str, paste(list.prefixes, collapse=" ")))



loadDat <- function(trait){

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
    # maybe keep the pvals? we seem to want these
    stan.obj <- extractDataStanMulti(list.ds2)
    return(stan.obj)
}




runM1 <- function(trait){
	# run model 1 for a specified trait

    dat <- loadDat(trait)
    print("Data loaded")
    dat$dat$K <- 2
    save(dat, file=sprintf("%s/dat_%s.RData", DATA.FOLDER, trait))
    print("LEARNING PARAMS")
    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    fit1 <- stan(file = "models/model1_no_loglik.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
    print("SAVING")
    rm(dat)
    print(fit1, pars=c("Sigma", "pi", "Omegacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    save(fit1, file=sprintf("%s/f_%s.RData", DATA.FOLDER, trait))

}

extractData <- function(trait){

    print("EXTRACTING")
    # load the data and fit to extract info about the run
    load(file=sprintf("%s/dat_%s.RData", DATA.FOLDER, trait))
    load(file=sprintf("%s/f_%s.RData", DATA.FOLDER, trait))
    
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
    write.table(next.row, file=sprintf("%s/summary_dat_%s_%s_%s.txt", DATA.FOLDER, trait, ndim, downsampled_str), row.names=FALSE, quote=FALSE)
}

calcErrBarsHerit <- function(trait){
    fit.file=sprintf("%s/f_%s.RData", DATA.FOLDER, trait)
    if (!file.exists(fit.file)){
            df <- data.frame(t(c("NA", trait, "NA", "NA")))
            colnames(df) <- c("value", "trait", "int", "sex")
            return(df)

        }
        load(file=fit.file)
        load(file=sprintf("%s/dat_%s.RData", DATA.FOLDER, trait))

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
	write.table(my.df3, file=sprintf("%s/h_err_%s_%s_%s.txt", DATA.FOLDER, trait, ndim, downsampled_str), row.names=FALSE, quote=FALSE, sep="\t")
    return(my.df3)
    }



   runM1(trait)
   extractData(trait)
   calcErrBarsHerit(trait)



