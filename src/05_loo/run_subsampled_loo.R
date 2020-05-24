require('rstan')
require('loo')
require('tictoc')
require('mnormt')

computeLLOne <- function(B, SE, p, sigmasq){
  zeros <- c(0,0)
  SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
  
  p_1 = p[1]*mnormt::dmnorm(B, zeros, SE_mat)
  p_2 = p[2]*mnormt::dmnorm(B, zeros, SE_mat + matrix(c(sigmasq[1], 0, 0, 0),2, 2))
  p_3 = p[3]*mnormt::dmnorm(B, zeros, SE_mat + matrix(c(0, 0, 0, sigmasq[2]),2,2))
  p_4 = p[4]*mnormt::dmnorm(B, zeros, SE_mat + matrix(c(sigmasq[3], 0, 0, sigmasq[4]), 2, 2))
  ll <- log(p_1)+log(p_2)+log(p_3)+log(p_4)
  return(ll)
}

computeLLOne_allDraws <- function(data_i, draws){
  B <- data_i[,1:2, drop=FALSE]
  SE <- data_i[,3:4, drop=FALSE]
  p <- draws[,1:4, drop=FALSE]
  s <- draws[,5:8, drop=FALSE]
  res <- lapply(1:nrow(draws), function(it) computeLLOne(B, SE, p[it,], s[it,]))
  return(unlist(res))
}

load("data/vary_priors_8/m2_fit_2.RData")
fit2 <- fit

ex_fit <- extract(fit2)
parameter_draws2 <- ex_fit$pi
parameter_draws2_ss <- ex_fit$sigmasq

load("data/vary_priors_8/dat_Testosterone.RData")
t_dat <- dat$dat

B.dat <- t_dat$B
SE.dat <- t_dat$SE


dat.df <- cbind(B.dat, SE.dat)
draws.df <- cbind(parameter_draws2, parameter_draws2_ss)
colnames(dat.df) <- c("b_f", "b_m", "se_f", "se_m")
colnames(draws.df) <- c("pi_1","pi_2","pi_3","pi_4",
                        "s_1", "s_2", "s_3", "s_4")

print(head(dat.df))
print(head(draws.df))

set.seed(4711)
tic()
loo_ss_2 <-
  loo_subsample(
    x=computeLLOne_allDraws,
    draws = draws.df,
    data = dat.df,
    observations = 200 # take a subsample of size 100
  )
print(toc())
print(loo_ss_2)
save(loo_ss_2, file="data/loo_ss_2.RData")