
args <- commandArgs(trailingOnly=TRUE)

model <- as.numeric(args[1])

trait <- args[2]
trait.type <- args[3] # binary or quant

if (length(args) >= 4){
   chrs=c(as.numeric(args[4]))
} else{
chrs <- c(1:22)
}

source('model_utils.R')
source('snp_utils.R')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


maf.cutoff <- 0.01 
se.cutoff <- 0.2


calcLoglik <- FALSE
filtDat <- TRUE

ndim <- 2
ndim_to_prefixes <- list("2"=c("zerosex", "onesex"), 
     "3"=c("pre_meno", "post_meno", "onesex"),
     "4"=c("under65_f", "over65_f", "under65_m", "over65_m"))

list.prefixes <- ndim_to_prefixes[[as.character(ndim)]]

require('tidyverse')

loadDat <- function(trait, trait.type){
    


    # for each trait
    list.ds <- lapply(list.prefixes, function(prefix) {
        all.dat <- do.call(rbind, lapply(chrs, function(chr) { getFile(prefix, chr, trait)}));
        colnames(all.dat)[1:3] <- c("CHR", "BP", "SNP");
        return(all.dat)
    })

    

    list.ds2 <- extractOverlappingRows(list.ds)
    
    if (filtDat){
      df.f <- list.ds2[[1]]
        df.m <- list.ds2[[2]]
	both.snps <- cbind(df.f[df.f$SNP  ,c("SNP", "CHR", "BP", "BETA","SE", "P")], 
             df.m[df.m$SNP ,c("BETA","SE", "P")])
       colnames(both.snps) <- c("SNP", "CHR", "BP", "B_f", "SE_f", "p_f", "B_m", "SE_m", "p_m")

       filt <- filter(both.snps,(p_f < 10**-2) | (p_m < 10**-2))
       list.filt <- list("1"=filter(df.f, SNP %in% filt$SNP), "2"=filter(df.m, SNP %in% filt$SNP))
       list.ds2 <- list.filt


    }


    stan.obj <- extractDataStanMulti(list.ds2)
    return(stan.obj)
}


runM2 <- function(trait, trait.type){
	# run model 2 for a specified trait

    dat <- loadDat(trait, trait.type)
    dat$dat$K <- 4
    save(dat, file=sprintf("%s/biomarker/m2/dat_%s.RData", DATA.FOLDER, trait))
    if (calcLoglik==TRUE){
        model.file <- "models/model2_loglik.stan" # this is v2...
    } else {
        model.file <- "models/model2.stan" # v2??
    }
    fit2 <- stan(file = model.file,  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/biomarker/m2/f_m2_%s.RData", DATA.FOLDER, trait))
}




runM3 <- function(trait, trait.type){
	# run model 2 for a specified trait

    dat <- loadDat(trait, trait.type)
    dat$dat$K <- 4
    save(dat, file=sprintf("%s/biomarker/m3/dat_%s.RData", DATA.FOLDER, trait))

    model.file <- "models/model3_v1.stan" # v2??

    fit2 <- stan(file = model.file,  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/biomarker/m3/f_m2_%s.RData", DATA.FOLDER, trait))
}


runM4 <- function(trait, trait.type){
	# run model 2 for a specified trait

    dat <- loadDat(trait, trait.type)
    dat$dat$K <- 6
    save(dat, file=sprintf("%s/biomarker/m4/dat_%s.RData", DATA.FOLDER, trait))

    model.file <- "models/model4_v1.stan" # v2??

    fit2 <- stan(file = model.file,  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/biomarker/m4/f_m2_%s.RData", DATA.FOLDER, trait))
}



runM2 <- function(trait, trait.type){
	# run model 2 for a specified trait

    dat <- loadDat(trait, trait.type)
    dat$dat$K <- 4
    save(dat, file=sprintf("%s/biomarker/m2/dat_%s.RData", DATA.FOLDER, trait))
    if (calcLoglik==TRUE){
        model.file <- "models/model2_loglik.stan" # this is v2...
    } else {
        model.file <- "models/model2.stan" # v2??
    }
    fit2 <- stan(file = model.file,  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
  
    print(fit2, pars=c("sigmasq", "pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/biomarker/m2/f_m2_%s.RData", DATA.FOLDER, trait))
}


runM2.a <- function(trait, trait.type){
	# run alternative model for model 2 
    
    dat <- loadDat(trait, trait.type)
    dat$dat$K <- 2

    if (calcLoglik==TRUE){
        model.file <- "models/model2_alt_loglik.stan"
    } else {
        model.file <- "models/model2_alt.stan"
    }

    fit2 <- stan(file = model.file,  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)

    print(fit2, pars=c("pi", "Sigma"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    print("SAVING")
    rm(dat)
    save(fit2, file=sprintf("%s/biomarker/m2/f_m2.a_%s.RData", DATA.FOLDER, trait))
}



extractSNPcat <- function(snp.df, df.f, df.m, category, trait, model_num=2){
    comp4 <- snp.df[which(snp.df$category==category),c("p1", "p2", "p3", "p4", "SNP")]
    if (length(comp4$SNP) > 0){
            both.snps <- cbind(df.f[df.f$SNP %in% comp4$SNP ,c("SNP", "CHR", "BP", "BETA","SE", "P")], 
             df.m[df.m$SNP %in% comp4$SNP,c("BETA","SE", "P")])
            colnames(both.snps) <- c("SNP", "CHR", "BP", "B_f", "SE_f", "p_f", "B_m", "SE_m", "p_m")
            both.snp.df <- both.snps[,c("SNP", "CHR", "BP", "B_f", "B_m", "SE_f", "SE_m", "p_m","p_f")] 
            both.snp.df <- merge(both.snp.df, comp4, by="SNP")       
    both.snp.df2 <- annotateSNP(both.snp.df)

    write.table(both.snp.df2, file=sprintf("%s/biomarker/m%s/snps%s_%s.txt", DATA.FOLDER, model_num, category, trait), row.names=FALSE)

    }   


}


### post-processing
extractData <- function(trait, model_num=2){
	print("Extracting")
    print(trait)

	load(file=sprintf("%s/biomarker/m%s/dat_%s.RData", DATA.FOLDER, model_num, trait)) 
    load(file=sprintf("%s/biomarker/m%s/f_m2_%s.RData", DATA.FOLDER, model_num, trait))

    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    Sigma <- getSigma(fit2)

    #write.table(data.frame(t(c(trait, unlist(p), unlist(sigmasq))), 
    #    file=sprintf("%s/biomarker/m%s/%s_summary.txt", DATA.FOLDER, model_num, trait), quote=FALSE, row.names=FALSE))

    # assign each SNP to a category
    posterior.df <- posteriorSNPtable(dat, fit2)
    write.table(posterior.df, file=sprintf("%s/biomarker/m%s/snp_table_%s.txt", DATA.FOLDER, model_num, trait), row.names=FALSE, quote=FALSE)
    print("Posterior table generated")

    # remove large files from the workspace
	rm(fit2)
    rm(dat)


    snp.df <- posterior.df

  cat.count <- table(snp.df$category)
    if ("2" %in% names(cat.count) | "3" %in% names(cat.count) | "4" %in% names(cat.count)){
        print(trait)
        print(cat.count)
        list.prefixes <- c("zerosex", "onesex")
        chrs <- c(1:22)
        list.ds <- lapply(list.prefixes, function(prefix) {
                all.dat <- do.call(rbind, lapply(chrs, function(chr) { getFile(prefix, chr, trait)}));
                colnames(all.dat)[1:3] <- c("CHR", "BP", "SNP");
                return(all.dat)
            })

        list.ds2 <- extractOverlappingRows(list.ds)
	df.f <- list.ds2[[1]]
	df.m <- list.ds2[[2]]

	extractSNPcat(snp.df, df.f, df.m, 2, trait, model_num)
	extractSNPcat(snp.df, df.f, df.m, 3, trait, model_num)
	extractSNPcat(snp.df, df.f, df.m, 4, trait, model_num)
    } else {
        print(sprintf("No sex-specific SNPs for %s", trait))
    }


}



if (model==2){
	runM2(trait, trait.type)
	extractData(trait)
} 
if (model == 1){
	runM2.a(trait, trait.type)
}

if (model == 3){
	runM3(trait, trait.type)
	extractData(trait, 3)
}

if (model == 4){
	runM4(trait, trait.type)
	extractData(trait, 4)

}
