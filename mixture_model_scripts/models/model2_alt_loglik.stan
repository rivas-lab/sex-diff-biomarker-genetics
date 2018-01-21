data {
    int<lower=1> K; // number of mixture components
    int<lower=1> N; // number of samples
    int<lower=1> M; // dimension, should be 2

    matrix[N, M] B;  // betas from observed data
    matrix[N, M] SE; // SE from observed data
}

transformed data{
    vector[M] zeros; //For now, means are zero
    matrix[M,M] SE_mat[N];
    zeros = rep_vector(0, M);

    // fill in the matrix of standard errors
    for (n in 1:N) {
        SE_mat[n] = diag_matrix(to_vector(SE[n]));
    }
}


parameters {
    simplex[K] pi; // mixing proportions, simplex????
    cholesky_factor_corr[M] L_Omega;
    vector<lower=0>[M] tau; 

}
transformed parameters{

    matrix[M,M] Sigma;
    matrix[M,M] Sigmas[K];

    Sigmas[1] = diag_matrix(rep_vector(0,2));
    Sigma = diag_pre_multiply(tau, L_Omega)*diag_pre_multiply(tau, L_Omega)';
    Sigmas[2] = Sigma;
}


model {
    vector[K] ps; // contributions of each
    // put in prior for sigmas
    tau ~ cauchy(0, 2.5);
    L_Omega ~ lkj_corr_cholesky(2);
    pi ~ dirichlet(rep_vector(1, K)); // may need to adjust
    
    for (n in 1:N){
        for (k in 1:K){
            ps[k]  = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, SE_mat[n] + Sigmas[k]);
        }
        target += log_sum_exp(ps);
    }
}

generated quantities {
    vector[K] ps;
    vector[N] log_lik;
    
    for (n in 1:N){
        for (k in 1:K){
           ps[k] = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, diag_matrix(to_vector(SE[n]))  + Sigma[k]);
        }
        log_lik[n] = log_sum_exp(ps);
    }
}
