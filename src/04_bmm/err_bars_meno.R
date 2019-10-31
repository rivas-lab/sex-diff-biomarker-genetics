
BASE.DIR <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/"

calcErrBarsHerit <- function(trait, indir, outdir, idx){
	
    require('rstan')
    source(sprintf("%s/mixture_model_scripts/model_utils.R", BASE.DIR))
    source(sprintf("%s/mixture_model_scripts/heritability_utils.R", BASE.DIR))
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

	indices <- c((16*(idx-1)+1):(16*(idx)))
	#indices <- c(1,2)
	res.df <- do.call(rbind, lapply(indices, function(i){ 
		p.i <- p[i,];
		s.i <- Sigma[i,,];
		dat.i <- labelCategories(dat, s.i, p.i);
        	h.i <- overallHeritability(dat.i, s.i, p.i);

		df <- data.frame(t(c(i, unlist(h.i) )))
		if (length(h.i)==3){
		   colnames(df) <- c("trait", "h_pre_post", "h_pre_m", "h_post_m")
		} else {
		   colnames(df) <- c("trait", "h_f", "h_m")
		}
	}))
	write.table(res.df, file=sprintf("%s/%s_%s.txt", outdir, trait, idx), row.names=FALSE, quote=FALSE)
    }

args = commandArgs(trailingOnly=TRUE)
trait <- args[1]
indir <- as.numeric(args[2])
outdir <- as.numeric(args[3])
idx <- as.numeric(args[4])
calcErrBarsHerit(trait, indir, outdir, idx)