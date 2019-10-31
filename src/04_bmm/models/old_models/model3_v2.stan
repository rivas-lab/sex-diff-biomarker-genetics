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
    vector<lower=0>[2] sigmasq;
    cholesky_factor_corr[M] L_Omega;
    cholesky_factor_corr[M] L_Theta;
    vector<lower=0>[M] tau; 

}
transformed parameters{

    matrix[M,M] Sigma;
    matrix[M,M] Sigmas[K];
    vector[2] a;
    vector[2] b;
    
    a[1] = sigmasq[1];
    a[2] = 0.0;
    b[1] = 0.0;
    b[2] = sigmasq[2];

    Sigmas[1] = diag_matrix(rep_vector(0,2));
    Sigmas[2] = diag_matrix(a);
    Sigmas[3] = diag_matrix(b);
    Sigma = diag_pre_multiply(tau, L_Omega)*diag_pre_multiply(tau, L_Omega)';
    Sigma = diag_pre_multiply(tau, L_Omega)*diag_pre_multiply(tau, L_Omega)';
    Sigmas[4] = Sigma;
    Sigmas[5] = Sigma2;
}


model {
    vector[K] ps; // contributions of each
    // put in prior for sigmas
    tau ~ cauchy(0, 2.5);
    L_Omega ~ lkj_corr_cholesky(2);
    L_Theta ~ lkj_corr_cholesky(2);
    sigmasq ~ inv_gamma(1,1);
    pi ~ dirichlet(rep_vector(1, K)); // may need to adjust
    
    for (n in 1:N){
        for (k in 1:K){
            ps[k]  = log(pi[k]) + multi_normal_cholesky_lpdf(B[n] | zeros, diag_pre_multiply(SE_mat[n], L_Theta) + Sigmas[k]);
        }
        target += log_sum_exp(ps);
    }
}

