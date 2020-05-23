# run_bmm.R
# E Flynn
# 11/17/2017
#
# Updated code for running BMM model for a trait.
#
# Key updates:
#   - sped up model/data loading by removing log_lik calculation because I was not using this
#   - extract heritability results as part of the analysis instead of post-processing
#   - extract results divided by chromosome (intended for X/XY/autosomal sub-analyses)

source('src/04_bmm/test_code/project_utils.R')

source('src/04_bmm/model_utils.R')
source('src/04_bmm/heritability_utils.R')
source('src/04_bmm/snp_utils.R')




require('R.utils')
require('data.table')
require('tidyverse')
require('reshape2')
require('parallel')


DATA.FOLDER <- "data/"
GWAS.DIR <- "data/gwas/"
snps.to.keep <- read.table(sprintf("%s/snp_filt_list_wX_v3.txt", DATA.FOLDER), header=FALSE, colClasses="character")

# PARSE ARGUMENTS
args = commandArgs(trailingOnly=TRUE) 

model <- as.numeric(args[1])
param_id <- args[2]
trait <- "Testosterone"
ndim <- 2
model_dir <- "tmp_models_sim2/"
out_dir <- "data/vary_priors_sim2/"
in_dir <- GWAS.DIR
suffix <- ""


se.cutoff <- 0.2

list.prefixes <- c("zerosex", "onesex")

# ---- list all the possible parameters ---- #
m2_params_list <- list(
    "1" =list("gamma"=1, "beta"=1),
    "2" =list("gamma"=1, "beta"=5),
    "3" =list("gamma"=0.1, "beta"=1),
    "4" =list("gamma"=10, "beta"=1),
    "5" =list("gamma"=1, "beta"=2),
    "6" =list("gamma"=0.1, "beta"=2),
    "7" =list("gamma"=0.1, "beta"=5)
)



m1_params_list <- list(
    "1"=list("tau"=2.5, "lkj"=2, "beta"=list(1,1)),
    "2"=list("tau"=2.5, "lkj"=2, "beta"=list(2,2)),
     "3"=list("tau"=2.5, "lkj"=2, "beta"=list(5,5)),
     "4"=list("tau"=2.5, "lkj"=2, "beta"=list(2,1)),
     "5"=list("tau"=2.5, "lkj"=2, "beta"=list(1,2)),
     "6"=list("tau"=2.5, "lkj"=0.5, "beta"=list(1,1)),
     "7"=list("tau"=2.5, "lkj"=1, "beta"=list(1,1)),
     "8"=list("tau"=2.5, "lkj"=1.5, "beta"=list(1,1)),
     "9"=list("tau"=2.5, "lkj"=2, "beta"=list(5,1)),
     "10"=list("tau"=2.5, "lkj"=2.5, "beta"=list(1,1)),
     "11"=list("tau"=2.5, "lkj"=3, "beta"=list(1,1)),
    "12"=list("tau"=0.5, "lkj"=2, "beta"=list(1,1)),
     "13"=list("tau"=25, "lkj"=2, "beta"=list(1,1))
)



write_m1_text <- function(params, param_id) {
    tau <- params$tau
    lkj <- params$lkj
    beta <- params$beta
    m1_text <- sprintf("data {
    int<lower=0> N; 
    int<lower=1> M; 
    matrix[N, M] B; 
    matrix[N, M] SE; 
    int<lower=1> K; 
}
transformed data{
    vector[M] zeros;
    matrix[M,M] SE_mat[N];
    zeros = rep_vector(0, M);

    for (n in 1:N) {
        SE_mat[n] = diag_matrix(to_vector(SE[n]));
    }
}

parameters {
    simplex[K] pi; 
    cholesky_factor_corr[M] L_Omega;
    vector<lower=0>[M] tau;
}

transformed parameters{
    matrix[M, M] Sigma;
    matrix[M, M] Sigmas[K];
    Sigma = diag_pre_multiply(tau, L_Omega)*diag_pre_multiply(tau, L_Omega)';

    Sigmas[1] = diag_matrix(rep_vector(0,M));
    Sigmas[2] = Sigma;
}

model {
    vector[K] ps; 
    pi ~ beta(%s, %s);
    tau ~ cauchy(0, %s);
    L_Omega ~ lkj_corr_cholesky(%s);


    for (n in 1:N){

       for (k in 1:K){
           ps[k] = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, SE_mat[n] + Sigmas[k]);
       }
       target += log_sum_exp(ps);

     }
}
generated quantities {
    vector[K] ps;
    vector[N] log_lik;
    matrix[M,M] Omegacor;

    Omegacor = multiply_lower_tri_self_transpose(L_Omega);
    
    for (n in 1:N){
        for (k in 1:K){
           ps[k] = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, SE_mat[n] + Sigmas[k]);
        }
        log_lik[n] = log_sum_exp(ps);
    }
}


",  beta[1], beta[2], tau, lkj)
    print(m1_text)
    sink(sprintf("%s/m1_%s.stan", model_dir, param_id))
    cat(m1_text)
    sink()
    }


write_m1_text2 <- function(params, param_id) {
    tau <- params$tau
    lkj <- params$lkj
    beta <- params$beta
    m1_text <- sprintf("data {
    int<lower=0> N; 
    int<lower=1> M; 
    matrix[N, M] B; 
    matrix[N, M] SE; 
    int<lower=1> K; 
}
transformed data{
    vector[M] zeros;
    matrix[M,M] SE_mat[N];
    zeros = rep_vector(0, M);

    for (n in 1:N) {
        SE_mat[n] = diag_matrix(to_vector(SE[n]));
    }
}

parameters {
    simplex[K] pi; 
    cholesky_factor_corr[M] L_Omega;
    vector<lower=0, upper=%s>[M] tau;
}

transformed parameters{
    matrix[M, M] Sigma;
    matrix[M, M] Sigmas[K];
    Sigma = diag_pre_multiply(tau, L_Omega)*diag_pre_multiply(tau, L_Omega)';

    Sigmas[1] = diag_matrix(rep_vector(0,M));
    Sigmas[2] = Sigma;
}

model {
    vector[K] ps; 
    pi ~ beta(%s,%s);
    L_Omega ~ lkj_corr_cholesky(%s);


    for (n in 1:N){

       for (k in 1:K){
           ps[k] = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, SE_mat[n] + Sigmas[k]);
       }
       target += log_sum_exp(ps);

     }
}
generated quantities {
    vector[K] ps;
    vector[N] log_lik;
    matrix[M,M] Omegacor;

    Omegacor = multiply_lower_tri_self_transpose(L_Omega);
    
    for (n in 1:N){
        for (k in 1:K){
           ps[k] = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, SE_mat[n] + Sigmas[k]);
        }
        log_lik[n] = log_sum_exp(ps);
    }
}

",  tau,beta[1], beta[2],lkj)
    print(m1_text)
    sink(sprintf("%s/m1_%s.stan", model_dir, param_id))
    cat(m1_text)
    sink()
    }

write_m2_text <- function(params, param_id) {
    gamma <- params$gamma
    beta <- params$beta

    m2_text <- sprintf("data {
    int<lower=1> K; 
    int<lower=1> N; 
    int<lower=1> M; 

    matrix[N, M] B; 
    matrix[N, M] SE;
}

transformed data{
    vector[M] zeros; 
    zeros = rep_vector(0, M);
}


parameters {
    simplex[K] pi; 
    vector<lower=0>[4] sigmasq;

}
transformed parameters{

    matrix[M,M] Sigma[K];
    vector[2] a;
    vector[2] b;
    vector[2] c;

    a[1] = sigmasq[1];
    a[2] = 0.0;
    b[1] = 0.0;
    b[2] = sigmasq[2];
    c[1] = sigmasq[3];
    c[2] = sigmasq[4];

    Sigma[1] = diag_matrix(rep_vector(0,2));
    Sigma[2] = diag_matrix(a);
    Sigma[3] = diag_matrix(b);
    Sigma[4] = diag_matrix(c);
}


model {
    vector[K] ps; 


    sigmasq ~ inv_gamma(%s,%s);
    pi ~ dirichlet(rep_vector(%s, K)); 

    for (n in 1:N){
        for (k in 1:K){
            ps[k]  = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, diag_matrix(to_vector(SE[n])) + Sigma[k]);
        }
        target += log_sum_exp(ps);
    }
}

",  gamma, gamma, beta)
    
    print(m2_text)
    
    sink(sprintf("%s/m2_%s.stan", model_dir, param_id))
    cat(m2_text, append=TRUE)
    sink()
}



write_m2_text2 <- function(params, param_id) {
    gamma <- params$gamma
    beta <- params$beta

    m2_text <- sprintf("data {
    int<lower=1> K; 
    int<lower=1> N; 
    int<lower=1> M; 

    matrix[N, M] B; 
    matrix[N, M] SE;
}

transformed data{
    vector[M] zeros; 
    zeros = rep_vector(0, M);
}


parameters {
    simplex[K] pi; 
    vector<lower=0>[4] sigmasq;
    vector[K] alpha;
    real<lower=0> kappa;
    simplex[K] theta[N];
}
transformed parameters{

    matrix[M,M] Sigma[K];
    vector[2] a;
    vector[2] b;
    vector[2] c;
    vector[K] alpha = kappa * phi;

    a[1] = sigmasq[1];
    a[2] = 0.0;
    b[1] = 0.0;
    b[2] = sigmasq[2];
    c[1] = sigmasq[3];
    c[2] = sigmasq[4];

    Sigma[1] = diag_matrix(rep_vector(0,2));
    Sigma[2] = diag_matrix(a);
    Sigma[3] = diag_matrix(b);
    Sigma[4] = diag_matrix(c);
}


model {
    vector[K] ps; 


    sigmasq ~ inv_gamma(%s,%s);
    pi ~ dirichlet(alpha); 

    for (n in 1:N){
        for (k in 1:K){
            ps[k]  = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, diag_matrix(to_vector(SE[n])) + Sigma[k]);
        }
        target += log_sum_exp(ps);
    }
}

",  gamma, gamma)
    
    print(m2_text)
    
    sink(sprintf("%s/m2_%s.stan", model_dir, param_id))
    cat(m2_text, append=TRUE)
    sink()
}




getData <- function(trait){
       # for each trait
    list.ds <- lapply(list.prefixes, function(prefix) {
        dat.1 <- fread(sprintf("%s/ukb24983_v2_hg19.%s_%s.genotyped.glm.linear", in_dir, trait, prefix), data.table=FALSE)
        colnames(dat.1) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P")
        colnames(dat.1)[1:3] <- c("CHR", "BP", "SNP");
			    
        rownames(dat.1) <- dat.1$SNP

        # remove NAs
        dat.2 <- dat.1[!is.na(dat.1$SE),]

        # SE filter
        dat.3 <- dat.2[dat.2$SE < QUANT.SE.CUTOFF,]

        # MAF filter
        dat.4 <- dat.3[dat.3$SNP %in% snps.to.keep$V1,]
        
        # downsample the data
        
        return(dat.4)

    })

    list.ds2 <- extractOverlappingRows(list.ds)

}

loadDat <- function(trait){
    N <- 10000 # number of samples
p <- c(0.8, 0.2) # fraction in each component (null, non-null)
S <- matrix(c(0.0012,0.0005, 0.0005, 0.0010),2,2)
fixS <- nearPD(S)$mat # nearest positive definite matrix
cov.data.M1.sim <- model1Sim(N, p, fixS)
    dat <- list("dat"=cov.data.M1.sim)
    return(dat)
}


loadDatM2 <- function(trait){
    N <- 10000 # number of samples
    p <- c(0.7, 0.08, 0.08, 0.14)
    sigmasq <- c(0.005, 0.001, 0.002, 0.004) 
    cov.data.M2.sim <- model2Sim(N, p, sigmasq)
    dat <- list("dat"=cov.data.M2.sim)
    return(dat)
}





runM1 <- function(trait, param_id){
    if (param_id == 14){
         write_m1_text2(m1_params_list[[1]], param_id)   
    } else {
         write_m1_text(m1_params_list[[param_id]], param_id)

    }
    

	# run model 1 for a specified trait
    dat <- loadDat(trait)
    dat$dat$K <- 2
    print(nrow(dat$dat$B))
    print("LEARNING PARAMS")
    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    set.seed(4444)    
    fit <- stan(file=sprintf("%s/m1_%s.stan", model_dir, param_id),  
            data = dat$dat,    
            chains = 4, warmup = 1000, iter = 2000, cores = 4, refresh = 200)
    rm(dat)
    print("saving")
    save(fit, file=sprintf("%s/m1_fit_%s.RData", out_dir, param_id))
    print(fit, pars=c("pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("extracting loo")
    require('loo')
    options(mc.cores = 4)
    log_lik1 <- extract_log_lik(fit)
    rm(fit)
    loo1 <- loo(log_lik1)
    print(loo1)
    save(loo1, file=sprintf("%s/loo_m1_%s.RData", out_dir, param_id))
    waic1 <- waic(log_lik1) 
    print(waic1)
    save(waic1, file=sprintf("%s/waic_m1_%s.RData", out_dir, param_id))
}


runM1a <- function(trait, param_id){
    write_m1_text(m1_params_list[[param_id]], param_id)

	# run model 1 for a specified trait
    dat <- loadDat(trait)
    dat$dat$K <- 2
    print(nrow(dat$dat$B))
    print("LEARNING PARAMS")
    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    set.seed(4444)
    fit <- stan(file="src/04_bmm/models/model1_loglik.stan",
            data = dat$dat,    
            chains = 4, warmup = 1000, iter = 2000, cores = 4, refresh = 200)
    rm(dat)
    print("saving")
    save(fit, file=sprintf("%s/m0_fit_%s.RData", out_dir, param_id))
    print(fit, pars=c("pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("extracting loo")
    require('loo')
    options(mc.cores = 4)
    log_lik1 <- extract_log_lik(fit)
    rm(fit)
    loo1 <- loo(log_lik1)
    print(loo1)
    save(loo1, file=sprintf("%s/loo_m1_%s.RData", out_dir, param_id))
    waic1 <- waic(log_lik1) 
    print(waic1)
    save(waic1, file=sprintf("%s/waic_m1_%s.RData", out_dir, param_id))
}


runM2 <- function(trait, param_id){
     if (param_id == 8){
         write_m2_text2(m2_params_list[[1]], param_id)   
    } else {
         write_m2_text(m2_params_list[[param_id]], param_id)

    }
    
    dat <- loadDatM2(trait)
    dat$dat$K <- 4

    
     print("LEARNING PARAMS")
    print(parallel::detectCores())
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    
    fit <- stan(file=sprintf("%s/m2_%s.stan", model_dir, param_id),  
            data = dat$dat,    
            chains = 4, warmup =1000, iter =2000, cores = 4, refresh = 200)
    rm(dat)
    print("saving")
    save(fit, file=sprintf("%s/m2_fit_%s.RData", out_dir, param_id))
    print("saved")
     require('loo')
    options(mc.cores = 4)
    log_lik1 <- extract_log_lik(fit)
    rm(fit)
    loo1 <- loo(log_lik1)
    print(loo1)
    save(loo1, file=sprintf("%s/loo_m2_%s.RData", out_dir, param_id))
    waic1 <- waic(log_lik1) 
    print(waic1)
    save(waic1, file=sprintf("%s/waic_m2_%s.RData", out_dir, param_id))
}

if (model==1){
   runM1(trait, param_id)

}
if (model==2){
    runM2(trait, param_id)
}

if (model==0){
  runM1a(trait, param_id)
}