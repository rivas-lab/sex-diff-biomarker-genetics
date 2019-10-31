# Model Utils
# Emily Flynn
# 10/10/2017
#
# Updated utilities for running prepping/loading data and running mixture models.


require('MASS')
require('Matrix')
require('mnormt')
require('qqman')
require('rstan')

# hard-coded cutoffs, we can adjust
BINARY.SE.CUTOFF <- 1 ## might want to adjust
QUANT.SE.CUTOFF <- 0.2

DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


# list of filtered variants - how was this constructed?? need to check


extractOverlappingRows <- function(list.ds){
    print(length(list.ds))
    print(dim(list.ds[[1]]))
    rows.to.keep <- rownames(list.ds[[1]])
    for (i in 2:length(list.ds)){
        dat <- list.ds[[i]]
        rows.to.keep <- intersect(rows.to.keep, rownames(dat))
    }
    list.ds2 <- lapply(list.ds, function(x) x[rows.to.keep,])

    return(list.ds2)
}

extractDataStanMulti <- function(all.dat){
    # put together betas and ses, sq se for SE matrix
    betas <- do.call(cbind, lapply(all.dat, function(x) x$BETA))
    ses <- do.call(cbind, lapply(all.dat, function(x) x$SE))
    pvals <- do.call(cbind, lapply(all.dat, function(x) x$P))
    se2 <- apply(ses, c(1,2), function(x) x^2)

    cov.data <- list(
        N = nrow(betas),
        M = length(all.dat),
        B = betas,
        SE = se2,
        K = 2
    )
    snps <- sapply(all.dat[[1]]$SNP, as.character)
    chr <- sapply(all.dat[[1]]$CHR, as.character)
    dat <- list("dat"=cov.data, "snp"=snps, "chr"=chr, "p"=pvals)
    return(dat)
}









#### ----- PARAMETER EXTRACTION ----- ####

### ----- MODEL 1 ----- ###

# extract Sigma, pi
getSigma <- function(fit){
    fit_summ_S <- summary(fit, pars=c("Sigma"), probs=c(0.05, 0.95))
    Sigma <- matrix(fit_summ_S$summary[,c("mean")], 2, 2)
    return(Sigma)
}

getPi <- function(fit){
    fit_summ_pi <- summary(fit, pars=c("pi"), probs=c(0.05, 0.95))
    p <- as.vector(fit_summ_pi$summary[,c("mean")])
    return(p)
}
getSigmaMulti <- function(fit1, ndim){
      fit_summ_S <- summary(fit1, pars=c("Sigma"), probs=c(0.05, 0.95))
      Sigma <- matrix(fit_summ_S$summary[,c("mean")], ndim, ndim)
      return(Sigma)
}


getRgMulti <- function(fit, ndim){
    #rg <- cov2cor(S)
    #return(rg[1,2])
    fit_summ_R <- summary(fit, pars=c("Omegacor"), probs=c(0.05, 0.50, 0.95))
    rg <- matrix(fit_summ_R$summary[,c("50%")], ndim, ndim)
    return(rg[upper.tri(rg, diag=FALSE)])
}


getRgConfMulti <- function(fit, ndim){ # 95% confidence interval
    fit_summ_R <- summary(fit, pars=c("Omegacor"), probs=c(0.025, 0.975))


    pairs <- combn(ndim, 2)
    list.labels <- apply(pairs, 2, function(x) sprintf("Omegacor[%s,%s]", x[1], x[2]))
    conf.l <- fit_summ_R$summary[list.labels,"2.5%"]
    conf.u <- fit_summ_R$summary[list.labels,"97.5%"]
    return(list("l"=conf.l, "u"=conf.u))
}




# get variances
getVars <- function(fit){
    fitS <- summary(fit, pars=c("sigmasq"), probs=c(0.05, 0.95))

    sigmasq <- as.vector(fitS$summary[,c("mean")])
    return(sigmasq)
}

