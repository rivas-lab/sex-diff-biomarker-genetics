# Model Utils
# Emily Flynn
# 10/10/2017
#
# Updated utilities for running prepping/loading data and running mixture models.


require('MASS')
require('Matrix')
require('mnormt')
require('qqman')
require('rstan')

# hard-coded cutoffs, we can adjust
BINARY.SE.CUTOFF <- 2 ## might want to adjust
QUANT.SE.CUTOFF <- 0.4

GWAS.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/results/"
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


# list of filtered variants
vars.to.keep <- read.table(sprintf("%s/snp_filt_list.txt", DATA.FOLDER), header=FALSE, 
    colClasses="character")


#### ----- DATA INPUT ----- ####

fileChecks <- function(file.f, file.m, chr, field){
	if (!file.exists(file.f)){
        print(sprintf("File missing for c%s trait:%s - female", chr, field))
        return(-1) # error
    } 
    if (!file.exists(file.m)){
        print(sprintf("File missing for c%s trait:%s - male", chr, field))
        return(-1)
    } 

    try.f <- try(read.delim(file.f))
    try.m <- try(read.delim(file.m))
    if (inherits(try.f, "try-error")){
        print(sprintf("Error loading files for c%s trait:%s - female", chr, field))
        return(-1)
    }
    if (inherits(try.m, "try-error")){
        print(sprintf("Error loading files for c%s trait:%s - male", chr, field))
        return(-1)
    }

   	return(1)
}

filtUkbDat <- function(d1, d2){
    
    # select only the rows with the additive model
    d1.1 <- d1[d1$TEST == "ADD",]
    d2.1 <- d2[d2$TEST == "ADD",]
    rownames(d1.1) <- d1.1$SNP
    rownames(d2.1) <- d2.1$SNP
    
    # select rows in both and reorder based on this
    joint.rows <- intersect(rownames(d1.1), rownames(d2.1))
    d1.2 <- d1.1[joint.rows,]
    d2.2 <- d2.1[joint.rows,]
    
    # remove NAs
    present.rows <- c((!is.na(d1.2$SE)) & (!is.na(d2.2$SE)))
    
    return(list('1'=d1.2[present.rows,], '2'=d2.2[present.rows,]))
}


getData <- function(chr, field){
    # load data for a chromosome and a data field. 
    prefix <- sprintf("%sukb24893_v2.%s", GWAS.FOLDER, field) 
    file.f <- paste(c(prefix, ".zerosex.PHENO1_c", chr, ".glm.linear.gz"), collapse="")
    file.m <- paste(c(prefix, ".onesex.PHENO1_c", chr, ".glm.linear.gz"), collapse="")

    my.classes = c("numeric", "numeric", "character", "factor", "factor", "character",
                   "numeric", "numeric", "numeric", "numeric", "numeric")

    checks <- fileChecks(file.f, file.m, chr, field)
    if (checks == -1){ return(NA) }
    
    dat.f <- read.delim(file.f, colClasses=my.classes)
    dat.m <- read.delim(file.m, colClasses=my.classes)
    filt.dat <- filtUkbDat(dat.f, dat.m)
    return(filt.dat)
}



getDataBin <- function(chr, field){
    # load data for a chromosome and a data field. 
    prefix <- sprintf("%sukb24893_v2.%s", GWAS.FOLDER, field) 
    file.f <- paste(c(prefix, ".zerosex.PHENO1_c", chr, ".glm.logistic.hybrid.gz"), collapse="")
    file.m <- paste(c(prefix, ".onesex.PHENO1_c", chr, ".glm.logistic.hybrid.gz"), collapse="")

    my.classes = c("numeric", "numeric", "character", "factor", "factor", "factor", "character",
                   "numeric", "numeric", "numeric", "numeric", "numeric") # this is diff from quant

    checks <- fileChecks(file.f, file.m, chr, field)
    if (checks == -1){ return(NA) }

    dat.f <- read.delim(file.f,  colClasses=my.classes)
    dat.m <- read.delim(file.m, colClasses=my.classes)

    filt.dat <- filtUkbDat(dat.f, dat.m)
    return(filt.dat)
}


filterMAF <- function(maf.cutoff){
    rem.snps <- read.table(sprintf("%s/snp_filt_metadata.txt", DATA.FOLDER), header=TRUE)
    filt.snps <- rem.snps[rem.snps$maf > maf.cutoff,]
    return(filt.snps$ID)
}

reformatData <- function(all.dat, trait.type, maf.cutoff=0.01){
	
	# check for missing files
    file.checks <- sapply(all.dat, function(x){
        ifelse (length(x) != 2, 0,
            ifelse( (ncol(x$`1`) > 2) & (ncol(x$`2`) > 2), 1, 0))
    })
    if (0 %in% file.checks){
        stop((paste(c("File loading issues for ", trait, ". Stopping now."), collapse="")))
	}

    # relabel some columns
    all.dat.f <- do.call(rbind, lapply(all.dat, function(x) x$`1`))
    colnames(all.dat.f)[1:3] <- c("CHR", "BP", "SNP")
    all.dat.m <- do.call(rbind, lapply(all.dat, function(x) x$`2`))
    colnames(all.dat.m)[1:3] <- c("CHR", "BP", "SNP")

    if (trait.type == "binary"){
	    all.dat.f$BETA <- log(all.dat.f$OR)
    	all.dat.m$BETA <- log(all.dat.m$OR)    	
    }

    # FILTER + REFORMAT
    if (maf.cutoff==0.01){
        snps.to.keep <- unlist(vars.to.keep[,1])

    } else { # if a non-default MAF is provided, filter by this
        snps.to.keep <- filterMAF(maf.cutoff)
    }
    dat.f <- all.dat.f[all.dat.f$SNP %in% snps.to.keep,]
    dat.m <- all.dat.m[all.dat.m$SNP %in% snps.to.keep,]        

    return(list('1'=dat.f, '2'=dat.m))
}



filterSE <- function(dat.f, dat.m, trait.type, cutoff='default'){
	# filter standard error

    if (cutoff=='default'){
        if (trait.type == "binary"){
            cutoff <- BINARY.SE.CUTOFF
        } 
        if (trait.type == "quant"){
            cutoff <- QUANT.SE.CUTOFF
        }        
    } # otherwise keep the provided cutoff


	low.se.rows <- c((dat.f$SE < cutoff) & (dat.m$SE < cutoff))
	filt.se.f <- dat.f[low.se.rows,]
	filt.se.m <- dat.m[low.se.rows,]
	
	return(list('1'=filt.se.f, '2'=filt.se.m))
}


extractDataStan <- function(filt.f, filt.m){
    # put together betas and ses, sq se for SE matrix
    betas <- cbind(filt.f$BETA, filt.m$BETA)
    ses <- cbind(filt.f$SE, filt.m$SE)
    se2 <- apply(ses, c(1,2), function(x) x^2)

    cov.data <- list(
        N = nrow(betas),
        M = 2,
        B = betas,
        SE = se2,
        K = 2
    )
    snps <- sapply(filt.f$SNP, as.character)
    dat <- list("dat"=cov.data, "snp"=snps)
    return(dat)
}

#### ----- MODEL RUNNING ----- ####


timeModel <- function(command){
    start.time <- Sys.time()
    res <- command
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    return(res)
}





#### ----- PARAMETER EXTRACTION ----- ####

### ----- MODEL 1 ----- ###

# extract Sigma, pi
getSigma <- function(fit){
    fit_summ_S <- summary(fit, pars=c("Sigma"), probs=c(0.05, 0.95))
    Sigma <- matrix(fit_summ_S$summary[,c("mean")], 2, 2)
    return(Sigma)
}

getPi <- function(fit){
    fit_summ_pi <- summary(fit, pars=c("pi"), probs=c(0.05, 0.95))
    p <- as.vector(fit_summ_pi$summary[,c("mean")])
    return(p)
}

# compute genetic correlation
getRg <- function(S){
    rg <- cov2cor(S)
    return(rg[1,2])
}




## Estimate SNP-level heritabilities 

getH <- function(B, SE, p, Sigma){
    zeros <-c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
    p_1 = p[1]*pmnorm(B, zeros, SE_mat)
    p_2 = p[2]*pmnorm(B, zeros, SE_mat + Sigma)
    prob_1 = p_1 / (p_1 + p_2)
    prob_2 = p_2 / (p_1 + p_2)
    category <- rbinom(1, 1, prob=prob_2) 
    if (category == 1){
        h_f <- Sigma[1,1] / (Sigma[1,1] + SE[1])
        h_m <- Sigma[2,2] / (Sigma[2,2] + SE[2])
        return(c(h_f, h_m))
    } else {
        return(c(NA,NA))
    }    
}

getPlotHeritabilities <- function(cov_dat, fit, type){
    B_dat <- cov_dat$B
    SE_dat <- cov_dat$SE
    N <- cov_dat$N
    Sigma <- getSigma(fit)
    p <- getPi(fit)
    hx <- sapply(1:N, function(i) getH(B_dat[i,], SE_dat[i,], p, Sigma))

    plot(hx[1,] ~ hx[2,], ylab="h_f", xlab="h_m", 
         main=paste("SNP-level heritabilities (", type, ")", sep=""), col="darkblue")
    lines(x=seq(0,1), y=seq(0,1), add=TRUE)
    return(hx)
}


# get variances
getVars <- function(fit){
    fitS <- summary(fit, pars=c("sigmasq"), probs=c(0.05, 0.95))

    sigmasq <- as.vector(fitS$summary[,c("mean")])
    return(sigmasq)
}

## Assign sample to a category based on estimated params

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




#### ----- PLOTTING ----- ####

mhPlot <- function(filt.f, filt.m, trait.name)  {
    png(filename=sprintf("%s/mh_plots_%s.png", DATA.FOLDER, trait.name), width=1200, height=1400)
    par(mfrow=c(2,1))
    ymax <- ceiling(-log(min(filt.f$P, filt.m$P), base=10))+1
    manhattan(filt.f, main=paste(trait.name, " (women)", sep=" "), ylim=c(0, ymax))
    manhattan(filt.m, main=paste(trait.name, " (men)", sep=" "), ylim=c(0, ymax))   
    dev.off()
}

BSEplot <- function(dat) { 
    plot(rbind(dat$B[,1], dat$B[,2]) ~ rbind(dat$SE[,1], dat$SE[,2]), main="", xlab="SE", ylab="BETA")
}

zplot <- function(dat, trait.name){
    plot(density(dat$B[,2]/dat$SE[,2]), col='blue', main = trait.name, xlab="B/SE")
    lines(density(dat$B[,1]/dat$SE[,1]), col='red')
}

