# modules
#require('rstan')
#require('qqman')
#require('loo')

DATA.DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data'
CODE.DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/mixture_model_scripts/'

args <- commandArgs(trailingOnly=TRUE)
trait <- args[1]

source(sprintf("%s/model_utils.R", CODE.DIR))

maf.range <- c(0.01, 0.05)
se.range <- c(0.4)

loadVaryPars <- function(trait){
    col.labels <- c("MAF", "SE", "N", "pi[1]", "pi[2]", "Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]", "rg")
    df <- data.frame(matrix(vector(),0, length(col.labels)), dimnames(c(0, col.labels)))

    for (i in maf.range){
        for (j in se.range){
            #print(sprintf("Looking at MAF: %s SE: %s", i, j))

            load(sprintf("%s/f_m1_%s_m%s_s%s.RData", DATA.FOLDER, trait, i, j))
                                        #print(dat$dat$N)
            m1.pi <- getPi(fit1)
            m1.Sigma <- getSigma(fit1)

            print(fit1, pars=c("Sigma", "pi", "Omegacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)                                        #print(m1.pi)
                                        #print(m1.Sigma)
            rg <- getRg(fit1)
                                        #print(rg)
            next.row <- c(i, j, dat$dat$N, unlist(m1.pi), unlist(m1.Sigma), rg)
            df <- rbind(df, next.row)
        }
    }

    colnames(df) <- col.labels
    return(df)
}

df <- loadVaryPars(trait)

write.table(df, file=sprintf('%s/%s_summary.txt', DATA.FOLDER, trait))

ggplot(data=full.df, aes(MAF, rg, colour = trait))+ geom_jitter(aes(shape=factor(SE)), size=2.5, width=0.002)+ scale_shape(solid = FALSE)
ggplot(data=full.df[full.df$SE==0.2,], aes(MAF, rg, colour = trait))+ xlim(0, 0.1)+geom_point()+ggtitle("Standard error cutoff 0.2")+ylim(0.68, 1.0)
ggplot(data=full.df[full.df$SE==0.4,], aes(MAF, rg, colour = trait))+ xlim(0, 0.1)+geom_point()+ggtitle("Standard error cutoff 0.4")+ylim(0.68, 1.0)

