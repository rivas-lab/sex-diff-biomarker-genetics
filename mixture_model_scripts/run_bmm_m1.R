# run_bmm_m1.R
# E Flynn
# 11/17/2017
#
# Updated code for running BMM model #1 a trait.
#
# Key updates:
#   - sped up model/data loading by removing log_lik calculation because I was not using this
#   - extract heritability results as part of the analysis instead of post-processing
#   - extract results divided by chromosome (intended for X/XY/autosomal sub-analyses)

source('model_utils.R')
source('heritability_utils.R')
require('R.utils')
require('data.table')

DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


# PARSE ARGUMENTS
args = commandArgs(trailingOnly=TRUE) # 1, 2, or 3 (2-alt)

model <- as.numeric(args[1])
#if (!model %in% c(1,2,3)){ stop ("please specify a model (1,2,or 3) for the first argument") }

trait <- args[2]
trait.type <- args[3] # binary or quant
if (!trait.type %in% c('binary', 'quant')){ stop ("please specify the type of trait (binary or quant") }


# fourth argument - number of dimensions
if (length(args > 3)){
    ndim <- as.numeric(args[4])
} else {
    ndim <- 2
}

# fifth argument - downsampled?
if (length(args) > 4){
    downsampled <- args[5]
} else { downsampled <- FALSE }

biomarker <- args[6]

maf.cutoff <- 0.01
se.cutoff <- 0.2

chrs <- c(1:22, "X", "XY") # make sure works for X, XY

ndim_to_prefixes <- list("2"=c("zerosex", "onesex"), 
    "3"=c("pre_meno", "post_meno", "onesex"),
    "4"=c("under65_f", "over65_f", "under65_m", "over65_m"))

list.prefixes <- ndim_to_prefixes[[as.character(ndim)]]

if (downsampled==TRUE){ 
    list.prefixes <- sapply(list.prefixes, function(prefix) paste(prefix, "d", sep="_"))
}

downsampled_str <- ifelse(downsampled, "downsampled", "")
print(sprintf("Running M1 for trait %s with %s dim, %s, with prefixes %s", trait, ndim, downsampled_str, paste(list.prefixes, collapse=" ")))

### generate a log file
#print(sprintf("Running M1 for trait %s with MAF cutoff and SE cutoff %s.", trait, maf.cutoff, se.cutoff))

#snps.to.keep <- filterMAF(maf.cutoff)


loadBiomarkerDat <- function(trait, sex){
  dir <- "/oak/stanford/groups/mrivas/projects/biomarkers/results/plink/combined"
  file.path <- sprintf("%s/%s_%s.glm.linear.gz", dir, trait, sex)
  dat <- fread(file.path, data.table=FALSE)

   # select only the rows with the additive model
    dat.1 <- dat[dat$TEST == "ADD",]
    

    # remove NAs
    dat.2 <- dat.1[!is.na(dat.1$SE),]

    # SE filter
    dat.3 <- dat.2[dat.2$SE < QUANT.SE.CUTOFF,]

    # MAF filter
    dat.4 <- dat.3[dat.3$ID %in% snps.to.keep$V1,]
    colnames(dat.4)[1:3] <- c("CHR", "BP", "SNP")
    rownames(dat.4) <- dat.4$SNP

    return(dat.4)
}

loadDat <- function(trait, trait.type){
    

    if (biomarker==TRUE){
     list.ds <- lapply(c("female", "male"), function(sex) loadBiomarkerDat(trait, sex))
    } else {
    # for each trait
    list.ds <- lapply(list.prefixes, function(prefix) {
        all.dat <- do.call(rbind, lapply(chrs, function(chr) { getFile(prefix, chr, trait)}));
        colnames(all.dat)[1:3] <- c("CHR", "BP", "SNP");
        return(all.dat)
    })
    }

    list.ds2 <- extractOverlappingRows(list.ds)

    stan.obj <- extractDataStanMulti(list.ds2)
    return(stan.obj)
}


runM1 <- function(trait, trait.type){
	# run model 1 for a specified trait

    dat <- loadDat(trait, trait.type)
    print("Data loaded")
    dat$dat$K <- 2
    save(dat, file=sprintf("%s/biomarker/dat_%s.RData", DATA.FOLDER, trait))
    print("LEARNING PARAMS")
    fit1 <- stan(file = "models/model1_no_loglik.stan",  
            data = dat$dat,    
            chains = 4, warmup = 200, iter = 600, cores = 4, refresh = 200)
    print("SAVING")
    rm(dat)
    print(fit1, pars=c("Sigma", "pi", "Omegacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    save(fit1, file=sprintf("%s/biomarker/f_%s.RData", DATA.FOLDER, trait))

}

extractData <- function(trait){

    print("EXTRACTING")
    # load the data and fit to extract info about the run
    load(file=sprintf("%s/biomarker/dat_%s.RData", DATA.FOLDER, trait))
    load(file=sprintf("%s/biomarker/f_%s.RData", DATA.FOLDER, trait))
    
    m1.pi <- getPi(fit1)
    m1.Sigma <- getSigmaMulti(fit1, ndim)
    rg <- getRgMulti(fit1, ndim)
    rg.c <- getRgConfMulti(fit1, ndim)

    rm(fit1)

    # assign each SNP to a component, estimate heritability
    dat <- labelCategories(dat, m1.Sigma, m1.pi) # label the SNPs with heritability
    h <- overallHeritability(dat, m1.Sigma, m1.pi)
    #hx <- getPlotHeritabilities(dat, m1.Sigma, trait)


    # estimate x, xy chromosome contributions to heritability
    #h.x <- getChrHeritability(dat, m1.Sigma, m1.pi, "X", trait)
    #h.xy <- getChrHeritability(dat, m1.Sigma, m1.pi, "XY", trait)
    #h.auto <- getChrHeritability(dat, m1.Sigma, m1.pi, "autosomal", trait)

    # write out the data
    next.row <- data.frame(t(c(trait, dat$dat$N, unlist(m1.pi), unlist(m1.Sigma), unlist(rg), unlist(rg.c$l), unlist(rg.c$u), unlist(h))))
        #unlist(h.x), unlist(h.xy), unlist(h.auto))))
    #col.labels <- c("trait", "N", "pi[1]", "pi[2]", "Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]", 
    #    "rg", "rg.l", "rg.u", "h") #, "h.x.f", "h.x.m", "h.xy.f", "h.xy.m", "h.auto.f", "h.auto.m")
    #colnames(next.row) <- col.labels

    # TODO - label these better

    write.table(next.row, file=sprintf("%s/biomarker/summary_dat_%s_%s_%s.txt", DATA.FOLDER, trait, ndim, downsampled_str), row.names=FALSE, quote=FALSE)
}

runM1(trait, trait.type)
extractData(trait)



