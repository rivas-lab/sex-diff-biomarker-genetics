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
    vector<lower=0>[M] sigmasame;
    positive_ordered[M] sigmaf;
    positive_ordered[M] sigmam;
    cholesky_factor_corr[M] L_Omega;
    vector<lower=0>[M] tau; 

}
transformed parameters{
    matrix[M,M] Sigma;
    matrix[M,M] Sigma2;
    matrix[M,M] Sigma3;
    matrix[M,M] Sigmas[K];

    vector[2] sigmafa;
    vector[2] sigmama;
    vector[2] sigmass;    

    sigmass[1] = sigmasame[1];
    sigmass[2] = sigmasame[1];
    sigmafa[1] = sigmaf[1];
    sigmafa[2] = sigmaf[2];
    sigmama[1] = sigmam[2];
    sigmama[2] = sigmam[1];
    Sigmas[1] = diag_matrix(rep_vector(0,2));

    Sigma = diag_pre_multiply(sigmass, L_Omega)*diag_pre_multiply(sigmass, L_Omega)';
    Sigmas[2] = Sigma;
    Sigma2 = diag_pre_multiply(sigmafa, L_Omega)*diag_pre_multiply(sigmafa, L_Omega)';
    Sigmas[3] = Sigma2;
    Sigma3 = diag_pre_multiply(sigmama, L_Omega)*diag_pre_multiply(sigmama, L_Omega)';
    Sigmas[4] = Sigma3;
}


model {
    vector[K] ps; // contributions of each
    // put in prior for sigmas
    tau ~ cauchy(0, 2.5);
    sigmasame ~ cauchy(0, 2.5);
    sigmaf ~ cauchy(0, 2.5);
    sigmam ~ cauchy(0, 2.5);
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
 	  matrix[M,M] Omegacor;
//    matrix[M,M] Thetacor;
    // matrix[N,K] p;
    Omegacor = multiply_lower_tri_self_transpose(L_Omega);
 //   Thetacor = multiply_lower_tri_self_transpose(L_Theta);

// for (n in 1:N){
// vector[K] p_raw;
// for (k in 1:K){
// p_raw[k] <- pi[k]*exp(multi_normal_cholesky_lpdf(B[n] | zeros, Sigmas[k] + diag_pre_multiply(SE_vec[n], L_Theta)));
// }
// for (k in 1:K){
// p[n,k] <- p_raw[k]/sum(p_raw);
// }
// }


   }
