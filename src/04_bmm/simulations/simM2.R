
require('rstan')
require('tidyverse')
source("model_utils.R")
source("snp_utils.R")
source("heritability_utils.R")
source("test_code/project_utils.R")



simRunM2 <- function(N, p, sigmasq, simID){
    
    cov.data.M2.sim <- model2Sim(N, p, sigmasq)
    fit2_sim <- stan(
      file = "models/model2.stan",  # Stan program
      data = cov.data.M2.sim,    # named list of data
      chains = 4,             # number of Markov chains
      warmup = 200,          # number of warmup iterations per chain
      iter = 600,            # total number of iterations per chain
      cores = 4,              
      refresh = 200          # show progress every 'refresh' iterations
      )
    p.est <- getPi(fit2_sim)
    sigmasq.est <- getVars(fit2_sim)

    if (length(sigmasq)==4) {
        Sigma.est <- matrix(c(sigmasq.est[3], 0, 0, sigmasq.est[4]),2,2)
    } else{
        Sigma.est <- getSigma(fit2_sim)
    }
 
    # calculate posteriors
    B_dat <- cov.data.M2.sim$B
    SE_dat <- cov.data.M2.sim$SE
    N <- cov.data.M2.sim$N
    
    posteriors <- lapply(1:N, function(i) computePosterior(B_dat[i,], SE_dat[i,], p.est, sigmasq.est, Sigma.est))
    posterior.df <- data.frame(do.call(rbind, posteriors))

    
    list.params <- list("p"=p, "s"=sigmasq, "p.est"=p.est, "sigmasq.est"=sigmasq.est)
    save(posterior.df, list.params, file=sprintf("../data/sim2_1028/m2sim_%s.RData", simID))
}

#pi.range <- seq(0.1, 0.9, 0.1)
samplePi <- function(){
    a <- runif(1, 0, 0.5)
    b <- runif(1, 0, (0.5-a)/2)
    c <- runif(1, 0, (0.5-a)/2)
    d <- 1-a-b-c
    return(c(d, b, c, a))
}



args <- commandArgs(trailingOnly=TRUE)
run_num <- args[1]

N <- 3000

sapply(1:10, function(idx){
	 p <- samplePi()
	 s.sel <- runif(3, 0.008,0.08)
	 s.sel <- c(s.sel[1],s.sel[2], s.sel[3], s.sel[3])
	 simID <- sprintf("%s_%s", run_num, idx)
	 print(simID)
	 simRunM2(N, p, s.sel, simID)
})