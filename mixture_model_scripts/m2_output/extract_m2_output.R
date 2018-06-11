
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"

source('model_utils.R')


extractData <- function(trait){

	print(trait)

	load(file=sprintf("%s/m2/dat_%s.RData", DATA.FOLDER, trait))
    load(file=sprintf("%s/m2/f_m2_%s.RData", DATA.FOLDER, trait))

    # fraction in non-null component
    p <- getPi(fit2)

    # sigmasq
    sigmasq <- getVars(fit2)
    print(fit2, pars=c("sigmasq", "pi"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
    write.table(data.frame(t(c(trait, unlist(p), unlist(sigmasq)))), file=sprintf("%s/m2/%s_summary.txt", DATA.FOLDER, trait), quote=FALSE, row.names=FALSE)
}

trait.list <- c('whr', 'trunk_fp', 'arm_fp', 'leg_fp', '49', '50', '4079')
sapply(trait.list, extractData)