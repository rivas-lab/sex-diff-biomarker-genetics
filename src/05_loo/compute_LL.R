require('rstan')
require('loo')
require('tictoc')
require('mnormt')


options(mc.cores = 4)


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


args <- commandArgs(trailingOnly=TRUE)
model_id <- args[1]
idx <- as.numeric(args[2])

load(sprintf("data/vary_priors_8/m2_fit_%s.RData", model_id))


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

#set.seed(2)
#my_s <- sample(1:nrow(dat.df), 20000)
#dat.df <- dat.df[my_s,]

# now let's index through
start.idx <- ((idx-1)*3000)+1
end.idx <- min(idx*3000, nrow(dat.df))

list.s <- start.idx:end.idx

tic()
ll_chunk <- lapply(list.s, function(i) computeLLOne_allDraws(dat.df[i,, drop=FALSE], draws.df))
toc()

#save(ll_chunk, file=sprintf("data/ll_chunk_%s.RData", idx))
save(ll_chunk, file=sprintf("data/full_ll/%s_ll_chunk_%s.RData", model_id, idx))