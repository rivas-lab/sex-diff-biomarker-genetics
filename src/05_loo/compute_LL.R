require('rstan')
require('loo')
require('tictoc')
require('mnormt')


options(mc.cores = 4)



computeLLM2 <- function(B, SE_sq, p, sigmasq){
  zeros <- c(0,0)
  SE_mat <- diag(SE_sq, nrow=2)
  
  p_1 = log(p[1])+mnormt::dmnorm(B, zeros, SE_mat, log=TRUE)
  p_2 = log(p[2])+mnormt::dmnorm(B, zeros, SE_mat + matrix(c(sigmasq[1], 0, 0, 0),2, 2), log=TRUE)
  p_3 = log(p[3])+mnormt::dmnorm(B, zeros, SE_mat + matrix(c(0, 0, 0, sigmasq[2]),2,2), log=TRUE)
  p_4 = log(p[4])+mnormt::dmnorm(B, zeros, SE_mat + matrix(c(sigmasq[3], 0, 0, sigmasq[4]), 2, 2), log=TRUE)
  ll <- log(exp(p_1)+exp(p_2)+exp(p_3)+exp(p_4))
  return(ll)
}

computeLLM1 <- function(B, SE_sq, p, Sigma){
  zeros <- c(0,0)
  SE_sq_mat <- diag(SE_sq, nrow=2) 

  p_1 = log(p[1])+mnormt::dmnorm(B, zeros, SE_sq_mat, log=TRUE)
  p_2 = log(p[2])+mnormt::dmnorm(B, zeros, SE_sq_mat+Sigma, log=TRUE)
  ll <- log(exp(p_1)+exp(p_2))
  return(ll)
}



computeLLM2_allDraws <- function(data_i, draws){
  B <- data_i[,1:2, drop=FALSE]
  SE <- data_i[,3:4, drop=TRUE]
  p <- draws[,1:4, drop=FALSE]
  s <- draws[,5:8, drop=FALSE]
  res <- lapply(1:nrow(draws), function(it) computeLLM2(B, SE, p[it,], s[it,]))
  return(unlist(res))
}

computeLLM1_allDraws <- function(data_i, draws){
  B <- data_i[,1:2, drop=FALSE]
  SE <- data_i[,3:4, drop=TRUE]
  p <- draws$p
  S <- draws$S
  res <- lapply(1:nrow(p), function(it) computeLLM1(B, SE, p[it,], S[it,,]))
  
  return(unlist(res))
}




args <- commandArgs(trailingOnly=TRUE)
model_id <- args[1]
param_id <- args[2]
idx <- as.numeric(args[3])

# load the data
load("data/vary_priors_freeze/dat_Testosterone.RData")
t_dat <- dat$dat
B.dat <- t_dat$B
SE.sq.dat <- t_dat$SE
dat.df <- cbind(B.dat, SE.sq.dat)
colnames(dat.df) <- c("b_f", "b_m", "se_f", "se_m")

# load the fit and then run the operations
load(sprintf("data/vary_priors_freeze/m%s_fit_%s.RData", model_id, param_id))
ex_fit <- extract(fit)
parameter_draws_p <- ex_fit$pi

# now let's index through
set.seed(2)
my_s <- sample(1:nrow(dat.df), 20000)
dat.df <- dat.df[my_s,]
list.s <- ((idx-1)*2000+1):(idx*2000)


if (model_id == 1){
  parameter_draws_S <- ex_fit$Sigma
  draws.df <- list("p"=parameter_draws_p, "S"=parameter_draws_S)
  tic()
  ll_chunk <- lapply(list.s, function(i) computeLLM1_allDraws(dat.df[i,, drop=FALSE], draws.df))
  print(toc())
  
}
if (model_id==2){
  parameter_draws_ss <- ex_fit$sigmasq
  draws.df <- cbind(parameter_draws_p, parameter_draws_ss)
  colnames(draws.df) <- c("pi_1","pi_2","pi_3","pi_4",
                          "s_1", "s_2", "s_3", "s_4")
  tic()
  ll_chunk <- lapply(list.s, function(i) computeLLM2_allDraws(dat.df[i,, drop=FALSE], draws.df))
  print(toc())
}

save(ll_chunk, file=sprintf("data/ll_attempt2/m%s_%s_ll_chunk_%s.RData", model_id, param_id, idx))





