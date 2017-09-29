## quant_model.R
## E Flynn
## 9/12/2017
##
## Code for running models 1 and 2 on a quantitative trait. 

args = commandArgs(trailingOnly=TRUE)
model <- as.numeric(args[1])

trait <- 21001


library(rstan)
require('qqman')
source('project_utils.R')
source('project_utils_add.R')

DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"

# LIST OF FILTERED VARIANTS
vars.to.keep <- read.table(sprintf("%s/snp_filt_list.txt", DATA.FOLDER), header=FALSE, 
	colClasses="character")

loadDat <- function(trait, plot=FALSE, se.filt=FALSE){

	# read in data, filter
	all.dat <- lapply(1:22, function(x){ getData(as.character(x), trait)})

	# relabel some columns
    all.dat.f <- do.call(rbind, lapply(all.dat, function(x) x$`1`))
    colnames(all.dat.f)[1:3] <- c("CHR", "BP", "SNP")
    all.dat.m <- do.call(rbind, lapply(all.dat, function(x) x$`2`))
    colnames(all.dat.m)[1:3] <- c("CHR", "BP", "SNP")

	# FILTER + REFORMAT
    filt.f <- all.dat.f[all.dat.f$SNP %in% unlist(vars.to.keep[,1]),]
    filt.m <- all.dat.m[all.dat.m$SNP %in% unlist(vars.to.keep[,1]),]
    
    # FILTER SE
    if (se.filt==TRUE){
        low.se.rows <- c((filt.f$SE < 0.4) & (filt.m$SE < 0.4))
        filt.f <- filt.f[low.se.rows,]
        filt.m <- filt.m[low.se.rows,]        
    }


    if (plot==TRUE){
    	mhPlot(filt.f, filt.m, trait)
    }

    # put together betas and ses, sq se for SE matrix
    betas <- cbind(filt.f$BETA, filt.m$BETA)
    ses <- cbind(filt.f$SE, filt.m$SE)
    se2 <- apply(ses, c(1,2), function(x) x^2)

    cov.data <- list(
        N = nrow(betas),
        M = 2,
        B = betas,
        SE = se2,
        K = 2
    )
    snps <- sapply(filt.f$SNP, as.character)
    dat <- list("dat"=cov.data, "snp"=snps)
    return(dat)
}

runM1 <- function(trait){

    dat <- loadDat(trait, se.filt=TRUE)

    start.time <- Sys.time()
    fit1 <- stan(
        file = "model1_loglik.stan",  
        data = dat$dat,    
        chains = 4, warmup = 200, iter = 300, cores = 4, refresh = 200)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    print(fit1, pars=c("Sigma", "pi"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    save(dat, fit1, file=paste(c(DATA.FOLDER, "f_m1_", trait, ".RData"), collapse=""))
}

runM2 <- function(trait){

    dat <- loadDat(trait, se.filt=TRUE)

    ### RUN M2
    dat$dat$K <- 4
    start.time <- Sys.time()
    fit2 <- stan(
        file = "model2_loglik.stan",  
        data = dat$dat,    
        chains = 4, warmup = 200, iter = 300, cores = 4, refresh = 200)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    print(fit2, pars=c("sigmasq", "pi"), probs=c(0.025, 0.5, 0.975), digits_summary=5)

    save(dat, fit2, file=paste(c(DATA.FOLDER, "f_m2_", trait, ".RData"), collapse=""))
}

runM2.a <- function(trait){
     dat <- loadDat(trait, se.filt=TRUE)

    ### RUN M2 alt
    start.time <- Sys.time()
    fit2 <- stan(
        file = "model2_alt_loglik.stan",  
        data = dat$dat,    
        chains = 4, warmup = 200, iter = 300, cores = 4, refresh = 200)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    print(fit2, pars=c("sigmasq", "pi"), probs=c(0.025, 0.5, 0.975), digits_summary=5)

    save(dat, fit2, file=paste(c(DATA.FOLDER, "f_m2.a_", trait, ".RData"), collapse=""))
}

mhPlot <- function(filt.f, filt.m, trait.name)	{
	png(filename=sprintf("%s/mh_plots_%s.png", DATA.FOLDER, trait.name), width=1200, height=1400)
	par(mfrow=c(2,1))
	ymax <- ceiling(-log(min(filt.f$P, filt.m$P), base=10))+1
	manhattan(filt.f, main=paste(trait.name, " (women)", sep=" "), ylim=c(0, ymax))
	manhattan(filt.m, main=paste(trait.name, " (men)", sep=" "), ylim=c(0, ymax))	
	dev.off()
}

BSEplot <- function(dat) { 
	plot(rbind(dat$B[,1], dat$B[,2]) ~ rbind(dat$SE[,1], dat$SE[,2]), main="", xlab="SE", ylab="BETA")
}

zplot <- function(dat, trait.name){
	plot(density(rbind(dat$B[,1], dat$B[,2])/rbind(dat$SE[,1], dat$SE[,2])), main = trait.name, xlab="B/SE")
	#plot(density(dat$B[,1]/dat$SE[,1]), main = sprintf("%s - F", trait.name))
	#plot(density(dat$B[,1]/dat$SE[,1]), main = sprintf("%s - M", trait.name))

}

if (model == 1){
    runM1(trait)
}
if (model==2){
    runM2(trait)
}
if (model==3){
    runM2.a(trait)
}


