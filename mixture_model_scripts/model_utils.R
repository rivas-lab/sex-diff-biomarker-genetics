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
BINARY.SE.CUTOFF <- 1 ## might want to adjust
QUANT.SE.CUTOFF <- 0.2

GWAS.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/age_sex_meno/" # TODO - update
#GWAS.FOLDER <- "/oak/stanford/groups/mrivas/projects/biomarkers/results/plink/combined" # for biomarker data
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


# list of filtered variants - how was this constructed?? need to check
snps.to.keep <- read.table(sprintf("%s/snp_filt_list.txt", DATA.FOLDER), header=FALSE, 
    colClasses="character")


#### ----- DATA INPUT ----- ####

fileChecks <- function(my.file, dat.source, chr, field){
    if (!file.exists(my.file)){
        print(sprintf("File missing for c%s trait:%s %s", chr, field, dat.source))
        return(-1) # error
    } 
   
    try.f <- try(read.table(my.file))
   
    if (inherits(try.f, "try-error")){
        print(sprintf("Error loading files for c%s trait:%s - %s", chr, field, dat.source))
        return(-1)
    }
    return(1)
}

### LOAD SOME OF THE data
# NOTE - this is only for quantitative
getFile <- function(dat.source, chr, field){
    prefix <- sprintf("%sukb24893_v2.%s", GWAS.FOLDER, field) 

    file.dat <- paste(c(prefix, ".", dat.source, ".PHENO1_c", chr, ".glm.linear.gz"), collapse="")

    my.classes = c("character", "numeric", "character", "character","character", "character",
                   "numeric", "numeric", "numeric", "numeric", "numeric")

    col.labels <- c("CHROM", "POS", "ID", "REF", "ALT1", "TEST", "OBS_CT", 
        "BETA", "SE", "T_STAT", "P")

    checks <- fileChecks(file.dat, dat.source, chr, field)
    if (checks == -1){ return(NA) }
    
    dat <- read.table(file.dat, colClasses=my.classes, header=FALSE)
    colnames(dat) <- col.labels

    # -- FILTERING -- #

    # select only the rows with the additive model
    dat.1 <- dat[dat$TEST == "ADD",]
    rownames(dat.1) <- dat.1$ID

    # remove NAs
    dat.2 <- dat.1[!is.na(dat.1$SE),]

    # SE filter
    dat.3 <- dat.2[dat.2$SE < QUANT.SE.CUTOFF,]

    # MAF filter
    dat.4 <- dat.3[dat.3$ID %in% snps.to.keep$V1,]

    return(dat.4)
}

extractOverlappingRows <- function(list.ds){
    print(length(list.ds))
    print(dim(list.ds[[1]]))
    rows.to.keep <- rownames(list.ds[[1]])
    for (i in 2:length(list.ds)){
        dat <- list.ds[[i]]
        rows.to.keep <- intersect(rows.to.keep, rownames(dat))
    }
    list.ds2 <- lapply(list.ds, function(x) x[rows.to.keep,])

    return(list.ds2)
}

extractDataStanMulti <- function(all.dat){
    # put together betas and ses, sq se for SE matrix
    betas <- do.call(cbind, lapply(all.dat, function(x) x$BETA))
    ses <- do.call(cbind, lapply(all.dat, function(x) x$SE))
    se2 <- apply(ses, c(1,2), function(x) x^2)

    cov.data <- list(
        N = nrow(betas),
        M = length(all.dat),
        B = betas,
        SE = se2,
        K = 2
    )
    snps <- sapply(all.dat[[1]]$SNP, as.character)
    chr <- sapply(all.dat[[1]]$CHR, as.character)
    dat <- list("dat"=cov.data, "snp"=snps, "chr"=chr)
    return(dat)
}


getSigmaMulti <- function(fit1, ndim){
      fit_summ_S <- summary(fit1, pars=c("Sigma"), probs=c(0.05, 0.95))
      Sigma <- matrix(fit_summ_S$summary[,c("mean")], ndim, ndim)
      return(Sigma)
}


getRgMulti <- function(fit, ndim){
    #rg <- cov2cor(S)
    #return(rg[1,2])
    fit_summ_R <- summary(fit, pars=c("Omegacor"), probs=c(0.05, 0.95))
    rg <- matrix(fit_summ_R$summary[,c("mean")], ndim, ndim)
    return(rg[upper.tri(rg, diag=FALSE)])
}


getRgConfMulti <- function(fit, ndim){ # 95% confidence interval
    fit_summ_R <- summary(fit, pars=c("Omegacor"), probs=c(0.025, 0.975))


    pairs <- combn(ndim, 2)
    list.labels <- apply(pairs, 2, function(x) sprintf("Omegacor[%s,%s]", x[1], x[2]))
    conf.l <- fit_summ_R$summary[list.labels,"2.5%"]
    conf.u <- fit_summ_R$summary[list.labels,"97.5%"]
    return(list("l"=conf.l, "u"=conf.u))
}


# fileChecks <- function(file.f, file.m, chr, field){
#	if (!file.exists(file.f)){
#         print(sprintf("File missing for c%s trait:%s - female", chr, field))
#         return(-1) # error
#     } 
#     if (!file.exists(file.m)){
#         print(sprintf("File missing for c%s trait:%s - male", chr, field))
#         return(-1)
#     } 

#     try.f <- try(read.table(file.f))
#     try.m <- try(read.table(file.m))
#     if (inherits(try.f, "try-error")){
#         print(sprintf("Error loading files for c%s trait:%s - female", chr, field))
#         return(-1)
#     }
#     if (inherits(try.m, "try-error")){
#         print(sprintf("Error loading files for c%s trait:%s - male", chr, field))
#         return(-1)
#     }

#	return(1)
# }

# filtUkbDat <- function(d1, d2){
    
    
#     # select rows in both and reorder based on this
#     joint.rows <- intersect(rownames(d1.1), rownames(d2.1))
#     d1.2 <- d1.1[joint.rows,]
#     d2.2 <- d2.1[joint.rows,]
    
    
#     return(list('1'=d1.2[present.rows,], '2'=d2.2[present.rows,]))
# }


# getDataQuant <- function(chr, field){
#     # load data for a chromosome and a data field,
#     prefix <- sprintf("%sukb24893_v2.%s", GWAS.FOLDER, field) 
#     file.f <- paste(c(prefix, ".zerosex.PHENO1_c", chr, ".glm.linear.gz"), collapse="")
#     file.m <- paste(c(prefix, ".onesex.PHENO1_c", chr, ".glm.linear.gz"), collapse="")

#     my.classes = c("character", "numeric", "character", "character","character", "character",
#                    "numeric", "numeric", "numeric", "numeric", "numeric")

#     col.labels <- c("CHROM", "POS", "ID", "REF", "ALT1", "TEST", "OBS_CT", 
#         "BETA", "SE", "T_STAT", "P")

#     checks <- fileChecks(file.f, file.m, chr, field)
#     if (checks == -1){ return(NA) }
    
#     dat.f <- read.table(file.f, colClasses=my.classes, header=FALSE)
#     dat.m <- read.table(file.m, colClasses=my.classes, header=FALSE)
#     colnames(dat.f) <- col.labels
#     colnames(dat.m) <- col.labels
#     filt.dat <- filtUkbDat(dat.f, dat.m)
#     return(filt.dat)
# }



# getDataBin <- function(chr, field){
#     # load data for a chromosome and a data field. 
#     prefix <- sprintf("%sukb24893_v2.%s", GWAS.FOLDER, field) 
#     file.f <- paste(c(prefix, ".zerosex.PHENO1_c", chr, ".glm.logistic.hybrid.gz"), collapse="")
#     file.m <- paste(c(prefix, ".onesex.PHENO1_c", chr, ".glm.logistic.hybrid.gz"), collapse="")

#     my.classes = c("character", "numeric", "character", "factor", "factor", "factor", "character",
#                    "numeric", "numeric", "numeric", "numeric", "numeric") # this is diff from quant

#     col.labels <- c("CHROM", "POS", "ID", "REF", "ALT", "FIRTH?", "TEST", "OBS_CT", 
#         "OR", "SE", "T_STAT", "P")

#     checks <- fileChecks(file.f, file.m, chr, field)
#     if (checks == -1){ return(NA) }

#     dat.f <- read.table(file.f,  colClasses=my.classes)
#     dat.m <- read.table(file.m, colClasses=my.classes)
#     colnames(dat.f) <- col.labels
#     colnames(dat.m) <- col.labels

#     filt.dat <- filtUkbDat(dat.f, dat.m)
#     return(filt.dat)
# }



# filterMAF <- function(maf.cutoff){
#     rem.snps <- read.table(sprintf("%s/snp_filt_metadata.txt", DATA.FOLDER), header=TRUE)
#     filt.snps <- sapply(rem.snps[rem.snps$maf >=maf.cutoff,]$ID, as.character)

#     # load other SNPs filtering data - X, XY, Y, MT
#     alt_chr <- read.table(sprintf("%s/chr_qc/alt_chr_qc_table.txt", DATA.FOLDER), header=TRUE)
#     filt.snps2 <- sapply(alt_chr[(alt_chr$keep==1 & alt_chr$MAF >=maf.cutoff),]$SNP, as.character)

#     filt.snps.full <- c(filt.snps, filt.snps2)
#     return(filt.snps.full)
# }

reformatData <- function(all.dat, trait.type, maf.cutoff=0.01){
	


    # deal with single chromosome - ex. just X
    if (length(all.dat)==1){
        all.dat.f <- data.frame(all.dat[[1]][[1]])
        all.dat.m <- data.frame(all.dat[[1]][[2]])
    } else {
        all.dat.f <- do.call(rbind, lapply(all.dat, function(x) x$`1`))
        all.dat.m <- do.call(rbind, lapply(all.dat, function(x) x$`2`))        
    }

    # relabel some columns
    colnames(all.dat.f)[1:3] <- c("CHR", "BP", "SNP")
    colnames(all.dat.m)[1:3] <- c("CHR", "BP", "SNP")

    if (trait.type == "binary"){
	    all.dat.f$BETA <- log(all.dat.f$OR)
    	all.dat.m$BETA <- log(all.dat.m$OR)    	
    }

    # filter by MAF
    snps.to.keep <- filterMAF(maf.cutoff)

    
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
    chr <- sapply(filt.f$CHR, as.character)
    dat <- list("dat"=cov.data, "snp"=snps, "chr"=chr)
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
getRg <- function(fit){
    #rg <- cov2cor(S)
    #return(rg[1,2])
    fit_summ_R <- summary(fit, pars=c("Omegacor"), probs=c(0.05, 0.95))
    rg <- matrix(fit_summ_R$summary[,c("mean")], 2, 2)
    return(rg[1,2])
}


getRgConf <- function(fit){ # 95% confidence interval
    fit_summ_R <- summary(fit, pars=c("Omegacor"), probs=c(0.025, 0.975))
    return(fit_summ_R$summary["Omegacor[1,2]",c("2.5%", "97.5%")])
}

# ---- MODEL 2 ---- #

getH <- function(B, SE, p, Sigma){
    zeros <-c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
    p_1 = p[1]*dmnorm(B, zeros, SE_mat)
    p_2 = p[2]*dmnorm(B, zeros, SE_mat + Sigma)
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

