source('model_utils.R')
source('heritability_utils.R')
require('stringr')
require('tidyr')
require('dplyr')
DATA.FOLDER <- "/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data/"

args <- commandArgs(trailingOnly=TRUE)
trait.idx <- as.numeric(args[1])

ref <- read.delim(sprintf("%s/multi_run_gwas_input.txt", DATA.FOLDER), header=FALSE)
ref2 <- read.delim(sprintf("%s/set_traits2.txt", DATA.FOLDER), header=FALSE)
ref2$trait.name <- gsub(" ", "_", str_trim(ref2$V3))
ref2$trait <- sapply(ref2$V2, as.character)
ref$trait.name <- gsub(" ", "_", str_trim(paste(ref$V4, ref$V5)))
ref$trait <- sapply(ref$V1, as.character)

ref.total <- rbind(ref[,c("trait", "trait.name")], ref2[,c("trait", "trait.name")])

NUM.BLOCK <- 8
if (trait.idx == 10){
   selected.traits <- ref.total$trait[((trait.idx-1)*NUM.BLOCK):nrow(ref.total)]
} else{
   selected.traits <- ref.total$trait[((trait.idx-1)*NUM.BLOCK):((trait.idx)*NUM.BLOCK)]

}

calcErrBarsHerit <- function(trait){
    fit.file=sprintf("%s/dat_set/f_%s.RData", DATA.FOLDER, trait)


    if (!file.exists(fit.file)){
        df <- data.frame(t(c("NA", trait, "NA", "NA")))
        colnames(df) <- c("value", "trait", "int", "sex")
        return(df)

    }
    load(file=fit.file)
    load(file=sprintf("%s/dat_set/dat_%s.RData", DATA.FOLDER, trait))
    
    # extract all estimate
    list_of_draws <- rstan::extract(fit1)
    pi.draws <- list_of_draws$pi
    p <- pi.draws
    s.draws <- list_of_draws$Sigma
    Sigma <- s.draws
    
    # extract lower + upper pi
    ordered.p <- p[order(p[,2]),] # ordering p by the non-null component 
    p.lower <- ordered.p[0.125*nrow(ordered.p),]
    p.upper <- ordered.p[0.975*nrow(ordered.p),]

    # extract lower + upper sigma
    ordered.S <- Sigma[order(Sigma[,1,1]),,]
    s.upper <- ordered.S[0.975*dim(Sigma)[1],,]
    s.lower <- ordered.S[0.125*dim(Sigma)[1],,]
    
    # recalculate SNP membership
    dat2 <- dat
    dat2$categories <- NULL
    dat.u <- labelCategories(dat2, s.upper, p.upper)
    dat.l <- labelCategories(dat2, s.lower, p.lower)

    h.up <- overallHeritability(dat.u, s.upper, p.upper)
    h.low <- overallHeritability(dat.l, s.lower, p.lower)
    res <- list("up"=h.up, "low"=h.low)
    
    # reformat into data frame
    my.df <- cbind(t(as.data.frame(res)), trait)

    my.df2 <- data.frame(cbind(my.df, rownames(my.df)))
    colnames(my.df2) <- c("value", "trait", "h.type")

    my.df2 <- my.df2 %>% separate( h.type, c("int", "sex"))
    rownames(my.df2) <- NULL
    return(my.df2)
}

list.tab <- lapply(selected.traits, calcErrBarsHerit)
full.tab <- do.call(rbind, list.tab)
write.table(full.tab, file=sprintf("%s/herit/errorBars_%s.txt", DATA.FOLDER, trait.idx), row.names=FALSE, col.names=FALSE)
