require('loo')

args <- commandArgs(trailingOnly=TRUE)
idx <- as.numeric(args[1])
my_f <- list.files(path="data/ll_chunks/",pattern=sprintf("%s_ll_*", idx))
list.loo <- lapply(my_f, function(f){
  load(sprintf("data/ll_chunks/%s", f))
  ll_df <- do.call(cbind, ll_chunk)
  return(ll_df)
})



ll_df2 <- do.call(cbind, list.loo)
r_eff <- relative_eff(ll_df2, chain_id=c(rep(1,1000),rep(2,1000), rep(3,1000),rep(4,1000)))
loo.obj <- loo(ll_df2, r_eff=r_eff)
save(loo.obj, file=sprintf("data/loo_fit/my_loo_m2_fit%s.RData", idx))