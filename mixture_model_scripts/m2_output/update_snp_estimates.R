

source('model_utils.R')
#source('heritability_utils.R')
source('snp_utils.R')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"


maf.cutoff <- 0.01 
se.cutoff <- 0.2
chrs <- c(1:22)

trait.type <- 'quant'
args <- commandArgs(trailingOnly=TRUE)
trait.idx <-  as.numeric(args[1])
traits <- c("trunk_fp", "arm_fp", "leg_fp", "whr", "4079", "49", "50", "20022", "30150")
trait <- traits[trait.idx]


extractData <- function(trait){
	print("Extracting")
    print(trait)

    load(file=sprintf("%s/m2_v4/dat_%s.RData", DATA.FOLDER, trait)) 
    load(file=sprintf("%s/m2_v4/f_m2_%s.RData", DATA.FOLDER, trait))

    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    Sigma <- getSigma(fit2)

    #write.table(data.frame(t(c(trait, unlist(p), unlist(sigmasq))), 
    #    file=sprintf("%s/m2_v2/%s_summary.txt", DATA.FOLDER, trait), quote=FALSE, row.names=FALSE))

    # assign each SNP to a category
    posterior.df <- posteriorSNPtable(dat, fit2)
    write.table(posterior.df, file=sprintf("%s/m2_v4/snp_table_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
    print("Posterior table generated")

    # remove large files from the workspace
    rm(fit2)
    rm(dat)

    # assumes quant
    all.dat <- lapply(chrs, function(x){ getDataQuant(as.character(x), trait)})

    # reformat data, remove rows that are not shared
    dat.reform <- reformatData(all.dat, trait.type, maf.cutoff)
    filt.f <- dat.reform$`1`
    filt.m <- dat.reform$`2`

    # filter by standard error
    dat.filt <- filterSE(filt.f, filt.m, trait.type, cutoff=se.cutoff)
    filt.f <- dat.filt$`1`
    filt.m <- dat.filt$`2`

    snp.tab <- sexSpecSNPtables(filt.f, filt.m, posterior.df) #sexSpecSNPtables(posterior.df$SNP, filt.f, filt.m, posterior.df$category)
    f.tab <- annotateSNP(snp.tab$'1')
    m.tab <- annotateSNP(snp.tab$'2')
    write.table(f.tab, file=sprintf("%s/m2_v4/f_spec_snp_tab_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
    write.table(m.tab, file=sprintf("%s/m2_v4/m_spec_snp_tab_%s.txt", DATA.FOLDER, trait), row.names=FALSE, quote=FALSE)
}

extractData(trait)