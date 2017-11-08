# test_m1_vary_pars.R
#
# E Flynn
# 10/16/2017
#
# Testing results for varying cutoffs for SE, MAF.

source('model_utils.R')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"

args <- commandArgs(trailingOnly=TRUE)
# read in the data
trait <- args[1]
print(trait)
trait.type <- 'quant'
all.dat <- lapply(1:22, function(x){ getData(as.character(x), trait)})

print("Data loaded")

testPars <- function(all.dat, trait, trait.type, maf.cut, se.cut){

	# vary the MAF cutoff
	#   - TODO: break up reformatData --> faster
	dat.reform <- reformatData(all.dat, trait.type, maf.cut) # reformat data, remove rows not shared, vary MAF
	filt.f <- dat.reform$`1`
	filt.m <- dat.reform$`2`

	print(sprintf("  Data reformatted, MAF filtered. %s rows", nrow(filt.m)))

	# filter by standard error
	# vary the SE cutoff
	dat.filt <- filterSE(filt.f, filt.m, trait.type, se.cut)
	filt2.f <- dat.filt$`1`
	filt2.m <- dat.filt$`2`
	print(sprintf("  SE filtered. %s rows.", nrow(filt2.f)))
	#print(nrow(filt2.m))

	dat <- extractDataStan(filt2.f, filt2.m)

	dat$dat$K <- 2

	m1 <- stan_model("models/model1_loglik.stan")
	f1 <- timeModel(vb(m1, dat$dat))#, hessian=TRUE))
	print(f1)

	save(dat, f1, file=sprintf("%s/test_pars/test_vb_%s_m%s_s%s.RData", DATA.FOLDER, trait, maf.cut, se.cut))
}

#maf.range <- c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10)
                                        #se.range <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4)
maf.range <- c(0.01, 0.05, 0.10)
se.range <- c(0.2)
for (i in maf.range){
	for (j in se.range){
		print(sprintf("Looking at MAF: %s SE: %s", i, j))
		testPars(all.dat, trait, trait.type, i, j)
	}
}
