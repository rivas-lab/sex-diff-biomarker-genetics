## test_bin_model.R
## E Flynn
## 9/12/2017
##
## Code for testing running models 1 and 2 on a binary trait. 


library(rstan)
require('qqman')
source('project_utils.R')
source('project_utils_add.R')

trait <- 'RH107'
OUTPUT.COLS <- c("trait", "frac", "pi1", "pi2", "Sigma11", "Sigma21", "Sigma12", "Sigma22", "rg", 
        "pi2.1", "pi2.2", "pi2.3", "pi2.4", "sigmasq1", "sigmasq2", "sigmasq3", "sigmasq4", "f.specific", "m.specific")

# LIST OF FILTERED VARIANTS
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data"

vars.to.keep <- read.table(sprintf("%s/snp_filt_list.txt", DATA.FOLDER), header=FALSE, 
    colClasses="character")

getDataBin <- function(chr, field){
    # load data for a chromosome and a data field. 
    prefix <- paste("/scratch/PI/mrivas/users/erflynn/sex_div_gwas/results/ukb24893_v2.", field, sep="")
    file.f <- paste(c(prefix, ".zerosex.PHENO1_c", chr, ".glm.logistic.hybrid.gz"), collapse="")
    file.m <- paste(c(prefix, ".onesex.PHENO1_c", chr, ".glm.logistic.hybrid.gz"), collapse="")

    my.classes = c("numeric", "numeric", "character", "factor", "factor", "factor", "character",
                   "numeric", "numeric", "numeric", "numeric", "numeric") # this is diff from quant

    if (!file.exists(file.f)){
        print(paste(c("File missing for c", chr, " ", field, " female"), collapse=""))
        return(NA)
    } 
    if (!file.exists(file.m)){
        print(paste(c("File missing for c", chr, " ", field, " male"), collapse=""))
        return(NA)
    } 



    try.f <- try(read.delim(file.f,  colClasses=my.classes))
    try.m <- try(read.delim(file.m, colClasses=my.classes))
    if (inherits(try.m, "try-error")){
        print(paste(c("Error loading files for c", chr, " HC", field, " male"), collapse=""))
        return(NA)
    }
    if (inherits(try.f, "try-error")){
        print(paste(c("Error loading files for c", chr, " HC", field, " female"), collapse=""))
        return(NA)
    }

    dat.f <- read.delim(file.f,  colClasses=my.classes)
    dat.m <- read.delim(file.m, colClasses=my.classes)

    filt.dat <- filtUkbDat(dat.f, dat.m)
    return(filt.dat)
}


filterSELogOdds <- function(all.dat){
    
    # relabel some columns
    all.dat.f <- do.call(rbind, lapply(all.dat, function(x) x$`1`))
    colnames(all.dat.f)[1:3] <- c("CHR", "BP", "SNP")
    all.dat.m <- do.call(rbind, lapply(all.dat, function(x) x$`2`))
    colnames(all.dat.m)[1:3] <- c("CHR", "BP", "SNP")

    # convert OR --> log ODDs for Betas
	all.dat.f$BETA <- log(all.dat.f$OR)
	all.dat.m$BETA <- log(all.dat.m$OR)

    # FILTER OUT STANDARD ERROR > 1
    low.se.rows <- c((all.dat.f$SE < 1) & (all.dat.m$SE < 1))
    filt.se.f <- all.dat.f[low.se.rows,]
    filt.se.m <- all.dat.m[low.se.rows,]

    return(list("1"=filt.se.f, "2"=filt.se.m))
}



all.dat <- lapply(1:22, function(x){ getDataBin(as.character(x), trait)})


    # check for missing files
    file.checks <- sapply(all.dat, function(x){
        ifelse (length(x) != 2, 0,
            ifelse( (ncol(x$`1`) == 12) & (ncol(x$`2`) == 12), 1, 0))
    })
    if (0 %in% file.checks){
        print(paste(c("File loading issues for ", trait, ". Stopping now."), collapse=""))
	}


    # relabel some columns
    all.dat.f <- do.call(rbind, lapply(all.dat, function(x) x$`1`))
    colnames(all.dat.f)[1:3] <- c("CHR", "BP", "SNP")
    all.dat.m <- do.call(rbind, lapply(all.dat, function(x) x$`2`))
    colnames(all.dat.m)[1:3] <- c("CHR", "BP", "SNP")

    all.dat.f$BETA <- log(all.dat.f$OR)
    all.dat.m$BETA <- log(all.dat.m$OR)

    # FILTER + REFORMAT
    filt.f <- all.dat.f[all.dat.f$SNP %in% unlist(vars.to.keep[,1]),]
    filt.m <- all.dat.m[all.dat.m$SNP %in% unlist(vars.to.keep[,1]),]
    
    # FILTER SE, plot


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

save(f1, f2, file="/home/erflynn/code/ukbb_mixture_models/test_bin_v2.RData")
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


save(f1.v, f2.v, file="/home/erflynn/code/ukbb_mixture_models/test_bin_v2_v.RData")