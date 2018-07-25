GWAS.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/results/"

BINARY.SE.CUTOFF <- 1 ## might want to adjust
QUANT.SE.CUTOFF <- 0.2



fileChecks2 <- function(my.file, dat.source, chr, field){
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

### LOAD SOME OF THE LEG_FP Data
getFile <- function(dat.source, chr, field){
    prefix <- sprintf("%sukb24893_v2.%s", GWAS.FOLDER, field) 

    file.dat <- paste(c(prefix, ".", dat.source, ".PHENO1_c", chr, ".glm.linear.gz"), collapse="")

    my.classes = c("character", "numeric", "character", "character","character", "character",
                   "numeric", "numeric", "numeric", "numeric", "numeric")

    col.labels <- c("CHROM", "POS", "ID", "REF", "ALT1", "TEST", "OBS_CT", 
        "BETA", "SE", "T_STAT", "P")

    checks <- fileChecks2(file.dat, dat.source, chr, field)
    if (checks == -1){ return(NA) }
    
    dat <- read.table(file.dat, colClasses=my.classes, header=FALSE)
    colnames(dat) <- col.labels
    return(dat)
}

filtUkbDat3 <- function(d1, d2, d3, trait.type, se.cutoff="default"){
    
    # select only the rows with the additive model
    d1.1 <- d1[d1$TEST == "ADD",]
    d2.1 <- d2[d2$TEST == "ADD",]
    d3.1 <- d3[d3$TEST == "ADD",]
    rownames(d1.1) <- d1.1$SNP
    rownames(d2.1) <- d2.1$SNP
    rownames(d3.1) <- d3.1$SNP
    
    # select rows in both and reorder based on this
    joint.rows <- intersect(intersect(rownames(d1.1), rownames(d2.1)), rownames(d3.1))
    d1.2 <- d1.1[joint.rows,]
    d2.2 <- d2.1[joint.rows,]
    d3.2 <- d3.1[joint.rows,]
    
    # remove NAs - filter SE
    if (se.cutoff=="default"){
        if (trait.type == "binary"){
            se.cutoff <- BINARY.SE.CUTOFF
        } 
        if (trait.type == "quant"){
            se.cutoff <- QUANT.SE.CUTOFF
        }        
    } 
    
    present.rows <- c((!is.na(d1.2$SE) & d1.2$SE < se.cutoff) & (!is.na(d2.2$SE) & d2.2$SE < se.cutoff) 
                      & (!is.na(d3.2$SE & d3.2$SE < se.cutoff)))
        
    return(list('1'=d1.2[present.rows,], '2'=d2.2[present.rows,], '3'=d3.2[present.rows,]))
}

getDataQuantMeno <- function(chr, field, trait.type, se.cutoff="default"){

    f1 <- getFile("pre_meno", chr, field)
    f2 <- getFile("post_meno", chr, field)
    f3 <- getFile("onesex", chr, field)

    filt.dat <- filtUkbDat3(f1, f2, f3, trait.type, se.cutoff)
    return(filt.dat)
}


reformatDat1 <- function(all.dat, idx, trait.type, maf.cutoff=0.01){
    # deal with single chromosome - ex. just X
    if (length(all.dat)==1){
        all.dat2 <- data.frame(all.dat[[1]][[idx]])
    } else {
        all.dat2 <- do.call(rbind, lapply(all.dat, function(x) x[[idx]]))
    }
    colnames(all.dat2)[1:3] <- c("CHR", "BP", "SNP")

    if (trait.type == "binary"){
        all.dat2$BETA <- log(all.dat2$OR)
    }

    # filter by MAF
    all.dat3 <- all.dat2[all.dat2$SNP %in% snps.to.keep,]
                                          
    # already filtered by SE --> return
    return(all.dat3)
}
                                          
reformDat3 <- function(all.dat, trait.type){
    (ndim <- length(all.dat[[1]]) )
    reform.dat <- lapply(1:ndim, function(idx) reformatDat1(all.dat, idx, trait.type))
    return(reform.dat)
}                                          

getSigma3 <- function(fit1, ndim){
	  fit_summ_S <- summary(fit1, pars=c("Sigma"), probs=c(0.05, 0.95))
	  Sigma <- matrix(fit_summ_S$summary[,c("mean")], ndim, ndim)
	  return(Sigma)
}

extractDataStan3 <- function(all.dat){
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
      