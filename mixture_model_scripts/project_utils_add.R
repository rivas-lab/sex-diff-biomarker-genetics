# E Flynn
# additional project utils



ld.filter <- function(ds){
    # LD-filter
    ds.f <- ds$`1`
    ds.m <- ds$`2`
    ld.data <- read.delim("/scratch/users/erflynn/ukbb_data/ld_pruning/all.prune.in", header=FALSE, 
                          colClasses="character")
    ds.f.ld <- ds.f[(ds.f$SNP %in% ld.data$V1),]
    ds.m.ld <- ds.m[(ds.m$SNP %in% ld.data$V1),]
    return(list("1"=ds.f.ld, "2"=ds.m.ld))
}



computePosterior <- function(B, SE, p, sigmasq){

    zeros <- c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)

    p_1 = p[1]*pmnorm(B, zeros, SE_mat)
    p_2 = p[2]*pmnorm(B, zeros, SE_mat + matrix(c(sigmasq[1], 0, 0, 0),2, 2))
    p_3 = p[3]*pmnorm(B, zeros, SE_mat + matrix(c(0, 0, 0, sigmasq[2]),2,2))
    p_4 = p[4]*pmnorm(B, zeros, SE_mat + matrix(c(sigmasq[3], 0, 0, sigmasq[4]), 2,2))
    p_tot = p_1 + p_2+ p_3 + p_4
    prob_1 = p_1 / p_tot
    prob_2 = p_2 / p_tot
    prob_3 = p_3 / p_tot
    prob_4 = p_4 / p_tot
    return(list(prob_1, prob_2, prob_3, prob_4))
}

getAllPosteriors <- function(cov.dat, fit){
	B.dat <- cov.dat$B
    SE.dat <- cov.dat$SE
    N <- cov.dat$N
    sigmasq <- getVars(fit)
    p <- getPi(fit) 
    
    posteriors <- lapply(1:N, function(i) computePosterior(B.dat[i,], SE.dat[i,], p, sigmasq))
    posterior.df <- data.frame(do.call(rbind, posteriors))
    return(posterior.df)
}

posteriorSNPtable <- function(dat, fit){
    # generate a posterior SNP table

	posterior.df <- getAllPosteriors(dat$dat, fit)

	# assign to the category with the maximum posterior
	posterior.df$category <- apply(posterior.df, 1, function(x){
		return(which.max(x))
	})
	posterior.df$SNP <- dat$snp
	colnames(posterior.df) <- c("p1", "p2", "p3", "p4", "category", "SNP")

	posterior.df <- data.frame(apply(posterior.df, c(1,2), unlist))
	return(posterior.df)
}


increasingPsample <- function(dat, cutoff.p=10**(-5), fracs, K=2){
	## sample with increasing fractions 
	## i.e. frac 0.10 contains all of 0.20 + more randomly sampled
	## used for testing randomizations

    dat.f <- dat$`1`
    dat.m <- dat$`2`
    topP.m <- dat.f[dat.f$P < cutoff.p,]
    topP.f <- dat.m[dat.m$P < cutoff.p,]
    top.snps <- union(topP.m$SNP, topP.f$SNP)
    print(length(top.snps))

    # top rows
    top.rows <- (dat.f$SNP %in% top.snps)
    top.f <- dat.f[top.rows,]
    top.m <- dat.m[top.rows,]

    # other rows
    other.f <- dat.f[-top.rows,]
    other.m <- dat.m[-top.rows,]  

    reorder.ind <- sample(1:nrow(other.f), nrow(other.f), replace=FALSE)
    
    list.out <- lapply(fracs, function(frac.top)
    {
	    # compute how many more 
	    num.rem <- ((1-frac.top)/frac.top) *length(top.snps)
	    
	    # randomly sample remaining percentage
	    other.ind <- reorder.ind[1:num.rem]
	    sampled.f <- rbind(top.f, other.f[other.ind,])
	    sampled.m <- rbind(top.m, other.m[other.ind,])
	    
	    # put together betas and ses, sq se for SE matrix
	    betas <- cbind(sampled.f$BETA, sampled.m$BETA)
	    ses <- cbind(sampled.f$SE, sampled.m$SE)
	    se2 <- apply(ses, c(1,2), function(x) x^2)

	    cov.data <- list(
	        N = nrow(betas),
	        M = 2,
	        B = betas,
	        SE = se2,
	        K = K
	    )
	    snps <- sapply(sampled.f$SNP, as.character)
	    return(list("dat"=cov.data, "snp"=snps))

    })
	return(list.out)   
}

###### --- EXTRA BELOW - not currently used --- ####
                                      
model2SimAlt <- function(N, p, sigmasq){
    # simulates data with no sex-specific components 
    #   - only samples from null + diag(sigmasq)
    # sigmasq are a vector of two variances
    
    zeros <- c(0,0)
    
    # sample squared SEs
    se2 <- simSE2(N)

    ### SAMPLE BETAS
    # M0
    n.m0 <- round(p[1]*N)
    se.m0 <- matrix(se2[1:(2*n.m0)], n.m0, 2)
    betas.m0 <- do.call(rbind, lapply(1:n.m0, function(x) mvrnorm(1, zeros, diag(se.m0[x,]))))


    # M1
    n.m1 <- N - n.m0
    se.m1 <- matrix(se2[(2*n.m0+1):(2*N)], n.m1, 2)
    betas.m1 <- do.call(rbind, lapply(1:n.m1, function(x) 
        mvrnorm(1, zeros, diag(se.m1[x,])+diag(c(sigmasq[1],sigmasq[2])))))


    # put together
    betas <- do.call(rbind, list(betas.m0, betas.m1))
    ses <- do.call(rbind, list(se.m0, se.m1))
        
    cov.data.k2.sim <- list(
        N = N,
        M = 2,
        B = betas,
        SE = ses,
        K = 2
    )
    return(cov.data.k2.sim)
}                                      
                           



computeLikelihood4 <- function(B, SE, p, sigmasq){
    zeros <- c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
 	
 	p_1 = p[1]*pmnorm(B, zeros, SE_mat)
    p_2 = p[2]*pmnorm(B, zeros, SE_mat + matrix(c(sigmasq[1], 0, 0, 0),2, 2))
    p_3 = p[3]*pmnorm(B, zeros, SE_mat + matrix(c(0, 0, 0, sigmasq[2]),2,2))
    p_4 = p[4]*pmnorm(B, zeros, SE_mat + matrix(c(sigmasq[3], 0, 0, sigmasq[4]), 2,2))
    p_vec = c(p_1, p_2, p_3, p_4)
    p_vec <- p_vec/sum(p_vec)
    p_vec = p_vec[p_vec !=0]
    l = sum(sapply(p_vec, log))
    return(l)
}


computeLikelihood2 <- function(B, SE, p, sigmasq){
    zeros <- c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
 	p_1 = p[1]*pmnorm(B, zeros, SE_mat)
    p_2 = p[2]*pmnorm(B, zeros, SE_mat + matrix(c(sigmasq[1], 0, 0, sigmasq[2]),2, 2))
    p_vec = c(p_1, p_2)
    p_vec <- p_vec/sum(p_vec)
    p_vec = p_vec[p_vec !=0]
    l = sum(sapply(p_vec, log))
    return(l)
}



bayesFactor2vs4 <- function(cov_dat, fit, fit_alt){
    B_dat <- cov_dat$B
    SE_dat <- cov_dat$SE
    N <- cov_dat$N
    sigmasq <- getVars(fit)
    p <- getPi(fit)
    sigmasq.a <- getVars(fit_alt)
    p.a <- getPi(fit_alt)

    l.m <- sum(sapply(1:N, function(i) computeLikelihood4(B_dat[i,], SE_dat[i,], p, sigmasq)))
    l.n <- sum(sapply(1:N, function(i) computeLikelihood2(B_dat[i,], SE_dat[i,], p.a, sigmasq.a)))
    print(l.m)
    print(l.n)
    return(l.m/l.n)
}

rgConf <- function(fit){
    # 95% confidence interval for the genetic correlation
    fit_summ_S <- summary(fit, pars=c("Sigma"), probs=c(0.025, 0.975))
    Sigma.l <- matrix(fit_summ_S$summary[,c("2.5%")], 2, 2)
    Sigma.u <- matrix(fit_summ_S$summary[,c("97.5%")], 2, 2)
    rg.l <- getRg(Sigma.l)
    rg.u <- getRg(Sigma.u)
    return(list("l"=min(rg.l, rg.u), "u"=max(rg.l, rg.u)))
}



 