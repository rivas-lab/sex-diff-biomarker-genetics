
require('rstan')
require('tidyverse')
require('MASS')
require('mnormt')

getPosterior2 <- function(B, SE, p, Sigma){
        # get the posterior probability for a SNP

    zeros <- rep(0, length(SE)) #c(0,0)
    SE_mat <- diag(SE) #matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
    p_1 = p[1]*dmnorm(B, zeros, SE_mat)
    p_2 = p[2]*dmnorm(B, zeros, SE_mat + Sigma)
    prob_1 = log(p_1) - log(p_1 + p_2)
    prob_2 = log(p_2) - log(p_1 + p_2)
    return(exp(prob_2))
}


getPosteriorVec <- function(trait){
    load(sprintf("../data/1019/m1/f_%s.RData", trait))
    load(sprintf("../data/1019/dat_%s.RData", trait))
    
    median_vals <- summary(fit1)$summary[,"50%"]
    p <- as.vector(median_vals[c("pi[1]", "pi[2]")])

    Sigma <- matrix(median_vals[c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]")], 2,2)
    B_dat <- dat$dat$B
    SE_dat <- dat$dat$SE
    N <- dat$dat$N

    posteriors <- sapply(1:N, function(i) getPosterior2(B_dat[i,], SE_dat[i,], p, Sigma))
    save(posteriors, file=sprintf("../data/tmp_posteriors2/post_%s.RData", trait))
    return(posteriors)
}

BIOMARKER.DIR <- "../data/1019/"
biomarker_traits <- list.files(BIOMARKER.DIR, pattern="*.RData")
biomarkers <- sapply(biomarker_traits, function(x) strsplit(strsplit(x, "dat_",fixed=TRUE)[[1]][[2]], ".RData")[[1]][[1]])
names(biomarkers) <- NULL

#getPosteriorVec("Testosterone")
post_vec_list <- lapply(biomarkers, getPosteriorVec)