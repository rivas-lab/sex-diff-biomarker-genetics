require('loo')
options(mc.cores=4)	
args <- commandArgs(trailingOnly=TRUE)

model <- as.numeric(args[1])
idx <- as.numeric(args[2])
my_f <- list.files(path="data/ll_attempt2/",pattern=sprintf("m%s_%s_ll_*", model, idx))
list.loo <- lapply(my_f, function(f){
  load(sprintf("data/ll_attempt2/%s", f))
  ll_df <- do.call(cbind, ll_chunk)
  return(ll_df)
})



ll_df2 <- do.call(cbind, list.loo)
r_eff <- relative_eff(ll_df2, chain_id=c(rep(1,1000),rep(2,1000), rep(3,1000),rep(4,1000)))
loo.obj <- loo(ll_df2, r_eff=r_eff)
save(loo.obj, file=sprintf("data/ll_comb/loo_m%s_fit%s.RData", model, idx))