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
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


# PARSE ARGUMENTS
args = commandArgs(trailingOnly=TRUE) # 1, 2, or 3 (2-alt)

model <- as.numeric(args[1])
#if (!model %in% c(1,2,3)){ stop ("please specify a model (1,2,or 3) for the first argument") }

trait <- args[2]
trait.type <- args[3] # binary or quant
if (!trait.type %in% c('binary', 'quant')){ stop ("please specify the type of trait (binary or quant") }


maf.cutoff <- 0.01
se.cutoff <- 0.2

### generate a log file
print(sprintf("Running M1 for trait %s with MAF cutoff and SE cutoff %s.", trait, maf.cutoff, se.cutoff))


loadDat <- function(trait, trait.type){
	# function for loading the data

	# load all the data
	if (trait.type == 'binary'){
		all.dat <- lapply(c(1:22, "X", "XY"), function(x){ getDataBin(as.character(x), trait)})
	} 
	if (trait.type == 'quant') {
		all.dat <- lapply(c(1:22, "X", "XY"), function(x){ getDataQuant(as.character(x), trait)})
	}

	# reformat data, remove rows that are not shared
    dat.reform <- reformatData(all.dat, trait.type, maf.cutoff)
    filt.f <- dat.reform$`1`
    filt.m <- dat.reform$`2`

    # filter by standard error
    dat.filt <- filterSE(filt.f, filt.m, trait.type, cutoff=se.cutoff)
    filt.f <- dat.filt$`1`
    filt.m <- dat.filt$`2`


    # extract dat in a format for stan input
    dat <- extractDataStan(filt.f, filt.m)

    return(dat)
}



runM1 <- function(trait, trait.type){
	# run model 1 for a specified trait
    
    dat <- loadDat(trait, trait.type)
    dat$dat$K <- 2
    save(dat, file=sprintf("%s/dat_set/dat_%s.RData", DATA.FOLDER, trait))
    print("LEARNING PARAMS")
    fit1 <- stan(file = "models/model1_no_loglik.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
    print("SAVING")
    rm(dat)
    print(fit1, pars=c("Sigma", "pi", "Omegacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    save(fit1, file=sprintf("%s/dat_set/f_%s.RData", DATA.FOLDER, trait))

}

extractData <- function(trait){

    print("EXTRACTING")
    # load the data and fit to extract info about the run
    load(file=sprintf("%s/dat_set/dat_%s.RData", DATA.FOLDER, trait))
    load(file=sprintf("%s/dat_set/f_%s.RData", DATA.FOLDER, trait))
    
    m1.pi <- getPi(fit1)
    m1.Sigma <- getSigma(fit1)
    rg <- getRg(fit1)
    rg.c <- getRgConf(fit1)
    rg.l <- rg.c[[1]]
    rg.u <- rg.c[[2]]

    rm(fit1)

    # assign each SNP to a component, estimate heritability
    dat <- labelCategories(dat, m1.Sigma, m1.pi) # label the SNPs with heritability
    h <- overallHeritability(dat, m1.Sigma, m1.pi)
    hx <- getPlotHeritabilities(dat, m1.Sigma, trait)


    # estimate x, xy chromosome contributions to heritability
    h.x <- getChrHeritability(dat, m1.Sigma, m1.pi, "X", trait)
    h.xy <- getChrHeritability(dat, m1.Sigma, m1.pi, "XY", trait)
    h.auto <- getChrHeritability(dat, m1.Sigma, m1.pi, "autosomal", trait)

    # write out the data
    next.row <- data.frame(t(c(trait, dat$dat$N, unlist(m1.pi), unlist(m1.Sigma), rg, rg.l, rg.u, unlist(h), unlist(h.x), unlist(h.xy), unlist(h.auto))))
    col.labels <- c("trait", "N", "pi[1]", "pi[2]", "Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]", 
        "rg", "rg.l", "rg.u", "h.f", "h.m", "h.x.f", "h.x.m", "h.xy.f", "h.xy.m", "h.auto.f", "h.auto.m")
    colnames(next.row) <- col.labels
    write.table(next.row, file=sprintf("%s/dat_set/summary_dat_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
}

runM1(trait, trait.type)
extractData(trait)



