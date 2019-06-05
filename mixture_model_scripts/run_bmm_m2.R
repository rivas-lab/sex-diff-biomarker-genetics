
args <- commandArgs(trailingOnly=TRUE)

model <- as.numeric(args[1])
#if (!model %in% c(1,2,3)){ stop ("please specify a model (1,2,or 3) for the first argument") }

trait <- args[2]
trait.type <- args[3] # binary or quant

source('model_utils.R')
#source('heritability_utils.R')
source('snp_utils.R')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


maf.cutoff <- 0.01 
se.cutoff <- 0.2
chrs <- c(1:22)

calcLoglik <- FALSE


ndim <- 2
ndim_to_prefixes <- list("2"=c("zerosex", "onesex"), 
     "3"=c("pre_meno", "post_meno", "onesex"),
     "4"=c("under65_f", "over65_f", "under65_m", "over65_m"))

list.prefixes <- ndim_to_prefixes[[as.character(ndim)]]

loadDat <- function(trait, trait.type){
    


    # for each trait
    list.ds <- lapply(list.prefixes, function(prefix) {
        all.dat <- do.call(rbind, lapply(chrs, function(chr) { getFile(prefix, chr, trait)}));
        colnames(all.dat)[1:3] <- c("CHR", "BP", "SNP");
        return(all.dat)
    })

    

    list.ds2 <- extractOverlappingRows(list.ds)

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



extractSNPcat <- function(snp.df, df.f, df.m, category, trait){
    comp4 <- snp.df[which(snp.df$category==category),c("p1", "p2", "p3", "p4", "SNP")]
    if (length(comp4$SNP) > 0){
            both.snps <- cbind(df.f[df.f$SNP %in% comp4$SNP ,c("SNP", "CHR", "BP", "BETA","SE", "P")], 
             df.m[df.m$SNP %in% comp4$SNP,c("BETA","SE", "P")])
            colnames(both.snps) <- c("SNP", "CHR", "BP", "B_f", "SE_f", "p_f", "B_m", "SE_m", "p_m")
            both.snp.df <- both.snps[,c("SNP", "CHR", "BP", "B_f", "B_m", "SE_f", "SE_m", "p_m","p_f")] 
            both.snp.df <- merge(both.snp.df, comp4, by="SNP")       
    both.snp.df2 <- annotateSNP(both.snp.df)

    write.table(both.snp.df2, file=sprintf("%s/biomarker/m2/snps%s_%s.txt", DATA.FOLDER, category, trait), row.names=FALSE)

    }   


}


### post-processing
extractData <- function(trait){
	print("Extracting")
    print(trait)

	load(file=sprintf("%s/biomarker/m2/dat_%s.RData", DATA.FOLDER, trait)) 
    load(file=sprintf("%s/biomarker/m2/f_m2_%s.RData", DATA.FOLDER, trait))

    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    Sigma <- getSigma(fit2)

    #write.table(data.frame(t(c(trait, unlist(p), unlist(sigmasq))), 
    #    file=sprintf("%s/biomarker/m2/%s_summary.txt", DATA.FOLDER, trait), quote=FALSE, row.names=FALSE))

    # assign each SNP to a category
    posterior.df <- posteriorSNPtable(dat, fit2)
    write.table(posterior.df, file=sprintf("%s/biomarker/m2/snp_table_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
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

	extractSNPcat(snp.df, df.f, df.m, 2, trait)
	extractSNPcat(snp.df, df.f, df.m, 3, trait)
	extractSNPcat(snp.df, df.f, df.m, 4, trait)
    } else {
        print(sprintf("No sex-specific SNPs for %s", trait))
    }


}



if (model==2){
	runM2(trait, trait.type)
	extractData(trait)
} else {
	runM2.a(trait, trait.type)
}


