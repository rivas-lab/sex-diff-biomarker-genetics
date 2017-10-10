## run_bmm.R
## E Flynn
## Updated: 10/10/2017
##
## Code for running bayesian mixture model on a binary or quantitative trait.
##   Relies on the code in model_utils for loading + running models on the data. 
##
##
## To run:
##    Rscript run_bmm.R <model=(1,2,3)> <trait> <trait-type=("binary" | "quant")> [<filter=(TRUE | FALSE)>] [<test-type=('opt' | 'vb')>]
##
## There are three possible models - model "3" is the alternate two-component model with no 
##    sex-specific components for M2, we use this for model comparison.
## Trait refers to phenotype ID, binary/quant is the type of trait.
## The optional arguments are filtering (default= TRUE), and which type of test model to run 
##    (optimizing or variational bayes).


source('model_utils.R')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


# PARSE ARGUMENTS
args = commandArgs(trailingOnly=TRUE) # 1, 2, or 3 (2-alt)

model <- as.numeric(args[1])
if (!model %in% c(1,2,3)){ stop ("please specify a model (1,2,or 3) for the first argument") }


trait <- args[2]
trait.type <- args[3] # binary or quant
if (!trait.type %in% c('binary', 'quant')){ stop ("please specify the type of trait (binary or quant") }


filt <- ifelse(length(args)>3 & args[4] == "FALSE", FALSE, TRUE) # whether to filter, optional (default TRUE)
test.type <- ifelse(length(args) == 5), args[5], 'NA') # parse test type ('vb' or 'opt'), optional


loadDat <- function(trait, trait.type, se.filt=FALSE){
	# function for loading the data

	# load all the data
	if (trait.type == 'binary'){
		all.dat <- lapply(1:22, function(x){ getDataBin(as.character(x), trait)})
	} 
	if (trait.type == 'quant') {
		all.dat <- lapply(1:22, function(x){ getData(as.character(x), trait)})
	}

	# reformat data, remove rows that are not shared
    dat.reform <- reformatData(all.dat, trait.type)
    filt.f <- dat.reform$`1`
    filt.m <- dat.reform$`2`

    # filter by standard error
    if (se.filt==TRUE){
        dat.filt <- filterSE(filt.f, filt.m, trait.type)
        filt.f <- dat.filt$`1`
        filt.m <- dat.filt$`2`
    }

    # extract dat in a format for stan input
    dat <- extractDataStan(filt.f, filt.m)

    return(dat)
}


runM1 <- function(trait, trait.type, filt=FALSE){
	# run model 1 for a specified trait

    dat <- loadDat(trait, trait.type, se.filt=filt)
    dat$dat$K <- 2
    fit1 <- timeModel(
        stan(file = "models/model1_loglik.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 300, cores = 4, refresh = 200))
   
    print(fit1, pars=c("Sigma", "pi"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    save(dat, fit1, file=paste(c(DATA.FOLDER, "f_m1_", trait, ".RData"), collapse=""))
}

runM2 <- function(trait, trait.type, filt=FALSE){
	# run model 2 for a specified trait

    dat <- loadDat(trait, trait.type, se.filt=filt)
    dat$dat$K <- 4

    fit2 <- timeModel(
        stan(file = "models/model2_loglik.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 300, cores = 4, refresh = 200))
  
    print(fit2, pars=c("sigmasq", "pi"), probs=c(0.025, 0.5, 0.975), digits_summary=5)

    save(dat, fit2, file=paste(c(DATA.FOLDER, "f_m2_", trait, ".RData"), collapse=""))
}

runM2.a <- function(trait, trait.type, filt=FALSE){
	# run alternative model for model 2 
    
    dat <- loadDat(trait, trait.type, se.filt=filt)
    dat$dat$K <- 2
    fit2 <- timeModel(
        stan(file = "models/model2_alt_loglik.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 300, cores = 4, refresh = 200))
 
    print(fit2, pars=c("sigmasq", "pi"), probs=c(0.025, 0.5, 0.975), digits_summary=5)

    save(dat, fit2, file=paste(c(DATA.FOLDER, "f_m2.a_", trait, ".RData"), collapse=""))
}


# run test models if specified
if (!isNA(test.type)){
	testModels(trait, trait.type, filt, test.type)
	exit()
}


# run models 
if (model==1){
    runM1(trait, trait.type, filt)
}
if (model==2){
    runM2(trait, trait.type, filt)
}
if (model==3){
    runM2.a(trait, trait.type, filt)
}
