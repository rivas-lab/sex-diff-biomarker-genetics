 
require('rstan')
require('tidyverse')
require('MASS')
require('mnormt')
source("snp_utils.R")

source("../mixture_model_scripts/snp_utils.R")
getPosteriorVec <- function(trait){
load(sprintf("../data/1020/m2/f_m2_%s.RData", trait))
   load(sprintf("../data/1020/dat_%s.RData", trait))

    median_vals <- summary(fit2)$summary[,"50%"]
    p <- median_vals[c("pi[1]", "pi[2]", "pi[3]", "pi[4]")]
    sigmasq <- median_vals[c("sigmasq[1]", "sigmasq[2]", "sigmasq[3]", "sigmasq[4]")]

    if (length(sigmasq)==4) {
        Sigma <- matrix(c(sigmasq[3], 0, 0, sigmasq[4]),2,2)
    } else{
        print("error")
    }
 
    # calculate posteriors
    B_dat <- dat$dat$B
    SE_dat <- dat$dat$SE
    N <- dat$dat$N
    
    posteriors <- lapply(1:N, function(i) computePosterior(B_dat[i,], SE_dat[i,], p, sigmasq, Sigma))

 posterior.df <- data.frame(do.call(rbind,posteriors))
colnames(posterior.df) <- c("p1", "p2", "p3", "p4")
      posterior.df$SNP <- dat$snp

        post.df <- posterior.df %>% filter(p2 >= 0.5 | p3 >= 0.5 | p4 >= 0.5) 

post.df2 <- data.frame(apply(post.df, c(1,2), function(x) unlist(x)))
  
#save(posterior_df, file=sprintf("../data/tmp_posteriors2/m2_post_%s.RData", trait))

p.df <- data.frame(dat$p )
    se.df <- data.frame(dat$dat$SE)
    beta.df <- data.frame(dat$dat$B)
    snp.df <- data.frame(dat$snp)
    chr.df <- data.frame(dat$chr)
    full.df <- do.call(cbind, list(snp.df, chr.df, beta.df, se.df, p.df))
    colnames(full.df) <- c("SNP", "CHR", "B.f", "B.m", "SE.f", "SE.m", "P.f", "P.m")
comb.df <- right_join(full.df,post.df2, by="SNP")
comb.df %>% write_csv(sprintf("../data/res_1024/m2/snp_table_%s.txt", trait))

}
args <- commandArgs(trailingOnly=TRUE)
trait <- args[1]
getPosteriorVec(trait)
