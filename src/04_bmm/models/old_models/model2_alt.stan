data {
    int<lower=1> K; // number of mixture components
    int<lower=1> N; // number of samples
    int<lower=1> M; // dimension, should be 2

    matrix[N, M] B;  // betas from observed data
    matrix[N, M] SE; // SE from observed data
}

transformed data{
    vector[M] zeros; //For now, means are zero
    zeros = rep_vector(0, M);
}


parameters {
    simplex[K] pi; // mixing proportions, simplex????
    vector<lower=0>[2] sigmasq;

}
transformed parameters{

    matrix[M,M] Sigma[K];

    vector[2] c;
    
    c[1] = sigmasq[1];
    c[2] = sigmasq[2];

    Sigma[1] = diag_matrix(rep_vector(0,2));
    Sigma[2] = diag_matrix(c);
}


model {
    vector[K] ps; // contributions of each
    // put in prior for sigmas
    sigmasq ~ inv_gamma(1,1);
    pi ~ dirichlet(rep_vector(1, K)); // may need to adjust
    
    for (n in 1:N){
        for (k in 1:K){
            ps[k]  = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, diag_matrix(to_vector(SE[n])) + Sigma[k]);
        }
        target += log_sum_exp(ps);
    }
}