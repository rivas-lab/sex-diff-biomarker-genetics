## test_quant_model.R
## E Flynn
## 9/12/2017
##
## Code for testing running models 1 and 2 on a quantitative trait. 

trait <- 21001


library(rstan)
require('qqman')
source('project_utils.R')
source('project_utils_add.R')

DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data"

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


dat <- loadDat(trait, plot=TRUE, se.filt=TRUE)
m1 <- stan_model("model1.stan")

start.time <- Sys.time()

f1 <- optimizing(m1, dat$dat, hessian=TRUE)
print(f1)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

dat$dat$K <- 4

m2 <- stan_model("model2.stan")

start.time <- Sys.time()

f2 <- optimizing(m2, dat$dat, hessian=TRUE)
print(f2)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

save(f1, f2, file="/home/erflynn/code/ukbb_mixture_models/test_quant_v2.RData")
dat$dat$K <- 2

start.time <- Sys.time()
f1.v <- vb(m1, dat$dat)
print(f1.v)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

dat$dat$K <- 4

start.time <- Sys.time()

f2.v <- vb(m2, dat$dat)
print(f2.v)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)


save(f1.v, f2.v, file="/home/erflynn/code/ukbb_mixture_models/test_quant_v2_v.RData")