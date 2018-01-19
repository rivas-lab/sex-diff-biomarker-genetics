
args <- commandArgs(trailingOnly=TRUE)

model <- as.numeric(args[1])
#if (!model %in% c(1,2,3)){ stop ("please specify a model (1,2,or 3) for the first argument") }

trait <- args[2]
trait.type <- args[3] # binary or quant

source('model_utils.R')
#source('heritability_utils.R')
source('snp_utils.R')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


maf.cutoff <- 0.01 
se.cutoff <- 0.2
chrs <- c(1:22)

loadDat <- function(trait, trait.type){
	# function for loading the data

	# load all the data
	if (trait.type == 'binary'){
		all.dat <- lapply(chrs, function(x){ getDataBin(as.character(x), trait)})
	} 
	if (trait.type == 'quant') {
		#all.dat <- lapply(c(1:22, "X", "XY"), function(x){ getDataQuant(as.character(x), trait)})
		all.dat <- lapply(chrs, function(x){ getDataQuant(as.character(x), trait)})

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



runM2 <- function(trait, trait.type){
	# run model 2 for a specified trait

    dat <- loadDat(trait, trait.type)
    dat$dat$K <- 4
    save(dat, file=sprintf("%s/m2_v3/dat_%s.RData", DATA.FOLDER, trait))
    fit2 <- stan(file = "models/model2_v2.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/m2_v3/f_m2_%s.RData", DATA.FOLDER, trait))
}

runM2.a <- function(trait, trait.type){
	# run alternative model for model 2 
    
    dat <- loadDat(trait, trait.type, maf.cutoff, filt, se.cutoff)
    dat$dat$K <- 2

    fit2 <- stan(file = "models/model2_alt.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)

    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/m2/f_m2.a_%s.RData", DATA.FOLDER, trait))
}


### post-processing
extractData <- function(trait){
	print("Extracting")

	load(file=sprintf("%s/m2_v3/dat_%s.RData", DATA.FOLDER, trait))
    load(file=sprintf("%s/m2_v3/f_m2_%s.RData", DATA.FOLDER, trait))

    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    Sigma <- getSigma(fit2)

    #write.table(data.frame(t(c(trait, unlist(p), unlist(sigmasq))), 
    #    file=sprintf("%s/m2_v2/%s_summary.txt", DATA.FOLDER, trait), quote=FALSE, row.names=FALSE))

    # assign each SNP to a category
    posterior.df <- posteriorSNPtable(dat, fit2)
    write.table(posterior.df, file=sprintf("%s/m2_v3/snp_table_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
    print("Posterior table generated")

    # remove large files from the workspace
	rm(fit2)
    rm(dat)

    # assumes quant
	all.dat <- lapply(chrs, function(x){ getDataQuant(as.character(x), trait)})

	# reformat data, remove rows that are not shared
    dat.reform <- reformatData(all.dat, trait.type, maf.cutoff)
    filt.f <- dat.reform$`1`
    filt.m <- dat.reform$`2`

    # filter by standard error
    dat.filt <- filterSE(filt.f, filt.m, trait.type, cutoff=se.cutoff)
    filt.f <- dat.filt$`1`
    filt.m <- dat.filt$`2`

    snp.tab <- sexSpecSNPtables(posterior.df$SNP, filt.f, filt.m, posterior.df$category)
    f.tab <- annotateSNP(snp.tab$'1')
    m.tab <- annotateSNP(snp.tab$'2')
	write.table(f.tab, file=sprintf("%s/m2_v3/f_spec_snp_tab_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
	write.table(m.tab, file=sprintf("%s/m2_v3/m_spec_snp_tab_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
}



if (model==2){
	runM2(trait, trait.type)
	extractData(trait)
} else {
	runM2.a(trait, trait.type)
}


