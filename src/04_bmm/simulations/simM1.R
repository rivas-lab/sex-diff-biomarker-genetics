
require('rstan')
require('tidyverse')
source("model_utils.R")
source("snp_utils.R")
source("heritability_utils.R")
#source("test_code/project_utils.R")

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



simSE2 <- function(N){
    
    #### SAMPLE STANDARD ERRORS from a half-normal distribution
    # parameters were selected to approximate true distribution of standard errors in the data
    se <- sapply(rnorm(N*5, 0.01, 0.02), abs)
    filt.se <- se[se <0.4 & se > 0.005]
    se2 <- sapply(filt.se, function(x) x^2)
    return(se2)
}


model1Sim <- function(N, p, Sigma){
    zeros <- c(0,0)

    # sample squared SEs
    se2 <- simSE2(N)

    ### SAMPLE BETAS FOR EACH MODEL
    # M0
    n.m0 <- floor(p[1]*N)
    se.m0 <- matrix(se2[1:(2*n.m0)], n.m0, 2)
    betas.m0 <- do.call(rbind, lapply(1:n.m0, function(x) mvrnorm(1, zeros, diag(se.m0[x,]))))

    # M1
    n.m1 <- N - n.m0
    se.m1 <- matrix(se2[(2*n.m0+1):(2*N)], n.m1, 2)
    betas.m1 <- do.call(rbind, lapply(1:n.m1, function(x) mvrnorm(1, zeros, diag(se.m1[x,])+Sigma)))

    # put together
    betas <- rbind(betas.m0, betas.m1)
    ses <- rbind(se.m0, se.m1)

    cov.data.sim <- list(
        N = N,
        M = 2,
        B = betas,
        SE = ses,
        K = 2
    )
    return(cov.data.sim)
}


simRun <- function(N, p_nn, S, simID){
    p <- c(1-p_nn, p_nn)
    
    fixS <- nearPD(S)$mat # nearest positive definite matrix
    cov.data.M1.sim <- model1Sim(N, p, fixS)
    fit1_sim <- stan(
      file = "models/model1_no_loglik.stan",  # Stan program
      data = cov.data.M1.sim,    # named list of data
      chains = 4,             # number of Markov chains
      warmup = 200,          # number of warmup iterations per chain
      iter = 600,            # total number of iterations per chain
      cores = 4,              
      refresh = 200          # show progress every 'refresh' iterations
      )
    p.est <- getPi(fit1_sim)
    s.est <- getSigma(fit1_sim)
    
    # calculate posteriors
    B_dat <- cov.data.M1.sim$B
    SE_dat <- cov.data.M1.sim$SE

    N <- cov.data.M1.sim$N
    posteriors <- sapply(1:N, function(i) getPosterior2(B_dat[i,], SE_dat[i,], p.est, s.est))
    
    list.params <- list("p"=p_nn, "S"=fixS, "p.est"=p.est[[2]], "S.est"=s.est)
    save(posteriors, list.params, file=sprintf("../data/tmp_posteriors/sim/sim%s.RData", simID))
}

#pi.range <- seq(0.1, 0.9, 0.1)
s.vals <- seq(0.0001,0.003, 0.0002) 

args <- commandArgs(trailingOnly=TRUE)
run_num <- args[1]
p_nn <- as.numeric(run_num)*0.1

N <- 1000

sapply(1:10, function(idx){
	 s.sel <- sample(s.vals, 3, replace=TRUE)
	 S <- matrix(c(s.sel[1], s.sel[2], s.sel[2], s.sel[3]), 2, 2)
	 simID <- sprintf("%s_%s", run_num, idx)
	 print(simID)
	 simRun(N, p_nn, S, simID)
})