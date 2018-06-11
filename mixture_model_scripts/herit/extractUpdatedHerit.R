

source('model_utils.R')
source('heritability_utils.R')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"

args = commandArgs(trailingOnly=TRUE)
trait <- args[1]


extractData <- function(trait){

    print("EXTRACTING")
    # load the data and fit to extract info about the run
    load(file=sprintf("%s/dat_set/dat_%s.RData", DATA.FOLDER, trait))
    load(file=sprintf("%s/dat_set/f_%s.RData", DATA.FOLDER, trait))
    
    m1.pi <- getPi(fit1)
    m1.Sigma <- getSigma(fit1)
    rg <- getRg(fit1)
    rg.c <- getRgConf(fit1)
    rg.l <- rg.c[[1]]
    rg.u <- rg.c[[2]]

    rm(fit1)

    # assign each SNP to a component, estimate heritability
    dat <- labelCategories(dat, m1.Sigma, m1.pi) # label the SNPs with heritability
    num1 <- table(dat$categories)[1]
    num2 <- length(dat$categories) - num1
    save(dat, file=sprintf("%s/dat_set/dat_v2_%s.RData", DATA.FOLDER, trait))

    h <- overallHeritability(dat, m1.Sigma, m1.pi)
    hx <- getPlotHeritabilities(dat, m1.Sigma, trait)


    # estimate x, xy chromosome contributions to heritability
    h.x <- getChrHeritability(dat, m1.Sigma, m1.pi, "X", trait)
    h.xy <- try({getChrHeritability(dat, m1.Sigma, m1.pi, "XY", trait)})
    if (class(h.xy)=="try-error"){
		h.xy <- c(0,0)
    }
    print(h.xy)
    #h.xy <- getChrHeritability(dat, m1.Sigma, m1.pi, "XY", trait)
    h.auto <- getChrHeritability(dat, m1.Sigma, m1.pi, "autosomal", trait)

    # write out the data
    next.row <- data.frame(t(c(trait, dat$dat$N, unlist(m1.pi), unlist(m1.Sigma), num1, num2, rg, rg.l, rg.u, unlist(h), unlist(h.x), unlist(h.xy), unlist(h.auto))))
    col.labels <- c("trait", "N", "pi[1]", "pi[2]", "Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]", 
        "num1", "num2", "rg", "rg.l", "rg.u", "h.f", "h.m", "h.x.f", "h.x.m", "h.xy.f", "h.xy.m", "h.auto.f", "h.auto.m")
    colnames(next.row) <- col.labels
    write.table(next.row, file=sprintf("%s/dat_set/summary_dat_%s_v2.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
}

extractData(trait)
