


computePosterior <- function(B, SE, p, sigmasq, Sigma){

    zeros <- c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)

    p_1 = p[1]*dmnorm(B, zeros, SE_mat)
    p_2 = p[2]*dmnorm(B, zeros, SE_mat + matrix(c(sigmasq[1], 0, 0, 0),2, 2))
    p_3 = p[3]*dmnorm(B, zeros, SE_mat + matrix(c(0, 0, 0, sigmasq[2]),2,2))
    p_4 = p[4]*dmnorm(B, zeros, SE_mat + Sigma)
    p_tot = p_1 + p_2+ p_3 + p_4
    prob_1 = exp(log(p_1) - log(p_tot))
    prob_2 = exp(log(p_2) - log(p_tot))
    prob_3 = exp(log(p_3) - log(p_tot))
    prob_4 = exp(log(p_4) - log(p_tot))
    return(list(prob_1, prob_2, prob_3, prob_4))
}

getAllPosteriors <- function(cov.dat, fit){
	B.dat <- cov.dat$B
    SE.dat <- cov.dat$SE
    N <- cov.dat$N
    sigmasq <- getVars(fit)
    p <- getPi(fit) 
    if (length(sigmasq)==4) {
        Sigma <- matrix(c(sigmasq[3], 0, 0, sigmasq[4]),2,2)
    } else{
        Sigma <- getSigma(fit)
    }
    
    posteriors <- lapply(1:N, function(i) computePosterior(B.dat[i,], SE.dat[i,], p, sigmasq, Sigma))
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




annotateSNP <- function(tab){
    DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"

    if (is.null(dim(tab))){return(tab)}
    gene.table <- read.table(sprintf("%s/snp_gene_table.txt", DATA.FOLDER), colClasses='character')
    combined.tab <- merge(tab, gene.table, by.x="SNP", by.y="snp", all.x=TRUE)
    return(combined.tab)
}  