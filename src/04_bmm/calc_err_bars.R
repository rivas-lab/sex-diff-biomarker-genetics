calcErrBarsHerit <- function(trait, indir, idx){
	
    require('rstan')
    source("src/04_bmm/model_utils.R")
    source("src/04_bmm/heritability_utils.R")
    require('tidyverse')
    ndim <- 2

    
    fit.file=sprintf("%s/m1/f_%s.RData", indir, trait)
        load(file=fit.file)
        load(file=sprintf("%s/dat_%s.RData", indir, trait))
	

        # extract all estimate
        list_of_draws <- rstan::extract(fit1)
        pi.draws <- list_of_draws$pi
        p <- pi.draws
        s.draws <- list_of_draws$Sigma
        Sigma <- s.draws

	indices <- c((40*(idx-1)+1):(40*(idx)))
	#indices <- c(1,2)
	res.df <- do.call(rbind, lapply(indices, function(i){ 
		p.i <- p[i,];
		s.i <- Sigma[i,,];
		dat.i <- labelCategories(dat, s.i, p.i);
        	h.i <- overallHeritability(dat.i, s.i, p.i);
		df <- data.frame("draw"=i, "h_f"=h.i[[1]], "h_m"=h.i[[2]])
	}))
	write.table(res.df, file=sprintf("data/h_err/%s_%s.txt", trait, idx), row.names=FALSE, quote=FALSE)
    }

args = commandArgs(trailingOnly=TRUE)
trait <- args[1]
data.dir <- as.numeric(args[2])
if (data.dir == 1){
    indir <- "data/"

} 
if (data.dir==2){
    indir <- "gwas/"

}


idx <- as.numeric(args[3])
calcErrBarsHerit(trait, indir, idx)