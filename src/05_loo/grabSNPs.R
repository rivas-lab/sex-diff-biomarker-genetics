
source('src/04_bmm/model_utils.R')
source('src/04_bmm/heritability_utils.R')
source('src/04_bmm/snp_utils.R')

require('R.utils')
require('data.table')
require('tidyverse')
require('reshape2')
require('parallel')


args <- commandArgs(trailingOnly=TRUE)
idx <- as.numeric(args[1])

load(sprintf("data/vary_priors_8/m2_fit_%s.RData", idx)); 
     fit2 <- fit
     rm(fit)
load("data/vary_priors_8/dat_Testosterone.RData")
    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    Sigma <- getSigma(fit2)
 posterior.df <- posteriorSNPtable(dat, fit2)
write.table(posterior.df, file=sprintf("data/snp_table_Testosterone_m2_%s.txt", idx), row.names=FALSE, quote=FALSE)


    cat.count <- table(posterior.df$category)
    print(cat.count)


p.df <- dat$p 
    se.df <- dat$dat$SE
    beta.df <- dat$dat$B
    snp.df <- data.frame(dat$snp, stringsAsFactors=FALSE)
    chr.df <- data.frame(dat$chr, stringsAsFactors=FALSE)
    full.df <- do.call(cbind, list(snp.df, chr.df, beta.df, se.df, p.df))
    colnames(full.df) <- c("SNP", "CHR", "B.f", "B.m", "SE.f", "SE.m", "P.f", "P.m")
    comb.df <- cbind(full.df, posterior.df %>% dplyr::select(p1, p2, p3, p4, category))
    non.null.snps <- comb.df %>% dplyr::filter(category %in% c(2,3,4))
write.table(non.null.snps, file=sprintf("data/sig_snp_table_Testosterone_m2_%s.txt", idx), row.names=FALSE, quote=FALSE)
