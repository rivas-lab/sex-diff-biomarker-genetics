# heritability_utils.R
# E Flynn
# 11/17/2017
#
# Utilities for extracting heritability information.

DATA.DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data'


getPosterior <- function(B, SE, p, Sigma){
	# get the posterior probability for a SNP

    zeros <- rep(0, length(SE)) #c(0,0)
    SE_mat <- diag(SE) #matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
    p_1 = p[1]*dmnorm(B, zeros, SE_mat)
    p_2 = p[2]*dmnorm(B, zeros, SE_mat + Sigma)
    prob_1 = log(p_1) - log(p_1 + p_2)
    prob_2 = log(p_2) - log(p_1 + p_2)
    return(exp(prob_2))
}

getCategory <- function(posterior){
	# get the category of a SNP given the posterior
	# we assign SNPs to the non-null component if the posterior probability is >= 0.7

    category <- ifelse(posterior >= 0.7, 2, 1) # is poster
    return(category)   
}

labelCategories <- function(dat, Sigma, p){
	# update the dat list to contain posterior probabilities and categories

    B_dat <- dat$dat$B
    SE_dat <- dat$dat$SE
    N <- dat$dat$N

	posteriors <- sapply(1:N, function(i) getPosterior(B_dat[i,], SE_dat[i,], p, Sigma))
	dat$posteriors <- posteriors
    categories <- sapply(posteriors, getCategory)
	dat$categories <- categories	
	return(dat)
}


getH <- function(SE, Sigma, category){
	#  get the individual heritability for a specific SNP

    if (category == 2){
        h_f <- Sigma[1,1] / (Sigma[1,1] + SE[1])
        h_m <- Sigma[2,2] / (Sigma[2,2] + SE[2])
        return(c(h_f, h_m))
    } else {
        return(c(NA,NA))
    }    
}

getPlotHeritabilities <- function(dat, Sigma, trait.name){
	# get and plot the individual SNP-level heritabilities

	N <- (dat$dat$N)
    hx <- sapply(1:N, function(i) getH(dat$dat$SE[i,], Sigma, dat$categories[i]))

    png(sprintf('%s/dat_set/%s_heritability.png', DATA.DIR, trait.name))
    plot(hx[1,] ~ hx[2,], ylab="h_f", xlab="h_m", 
         main=paste("SNP-level heritabilities (", trait.name, ")", sep=""), col="darkblue")
    lines(x=seq(0,1), y=seq(0,1))
    dev.off()
    return(hx)
}

calcHeritability <- function(se.p2, Sigma, p){
	n <- nrow(se.p2)
    num_f <- n*(p[2])*Sigma[1,1]
    num_m <- n*(p[2])*Sigma[2,2]
    
    h_f <- num_f/(num_f + sum(se.p2[,1]))
    h_m <- num_m/(num_m + sum(se.p2[,2]))
    return(list("1"=h_f, "2"=h_m))

}

calcHeritabilityMulti <- function(se.p2, Sigma, p, idx){
    n <- nrow(se.p2)
    num_i <- n*(p[2])*Sigma[idx,idx]
    
    h_i <- num_i/(num_i + sum(se.p2[,idx]))
    return(h_i)
}

overallHeritability <- function(dat, Sigma, p){
	# compute the overall heritability for a quantitative trait

    se.p2 <- dat$dat$SE[dat$categories==2,]

    list.h <- sapply(1:nrow(Sigma), function(idx) calcHeritabilityMulti(se.p2, Sigma, p, idx))
    return(list.h)
}

### TODO - x chromosome specific

getChrHeritability <- function(dat, Sigma, p, chr, trait.name){

	if (chr == "autosomal"){
		chr.list = sapply(seq(1,22), as.character)
	} else {
		chr.list = c(chr)
	}

	chr.snps <- (dat$chr %in% chr.list)
	chr.categories <- dat$categories[chr.snps]
	se.chr <- dat$dat$SE[chr.snps,] 
	se.p2.chr <- se.chr[chr.categories==2,]
	num.snps <- nrow(se.chr)
    n <- nrow(se.p2.chr)

	print(sprintf("Out of %s %s chromosome SNPs, %s were assigned to the non-null component (%s).", num.snps, chr, n, n/num.snps))

    # get/plot SNP specific heritability

    # create a new data object with enough info
    dat.dat.new <- list(SE=se.chr, N=num.snps)
    dat.filt <- list(dat=dat.dat.new, categories=chr.categories)
    hx <- getPlotHeritabilities(dat.filt, Sigma, sprintf("%s_%s", trait.name, chr) )

	# calculate the overall heritability
    return(calcHeritability(se.p2.chr, Sigma, p))
}