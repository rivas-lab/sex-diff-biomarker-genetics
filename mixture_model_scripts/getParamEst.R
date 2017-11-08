# getParamEst.R
# E Flynn
# 10/24/2017
#
# Extract results from Model 1 

maf <- 0.01
se <- 0.2


DATA.DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/data'
CODE.DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/mixture_model_scripts/'

source(sprintf("%s/model_utils.R", CODE.DIR))



getH <- function(B, SE, p, Sigma){
    zeros <-c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
    p_1 = p[1]*pmnorm(B, zeros, SE_mat)
    p_2 = p[2]*pmnorm(B, zeros, SE_mat + Sigma)
    prob_1 = p_1 / (p_1 + p_2)
    prob_2 = p_2 / (p_1 + p_2)
    category <- ifelse(prob_2 >= 0.8, 2, 1) #rbinom(1, 1, prob=prob_2) 
    if (category == 2){
        h_f <- Sigma[1,1] / (Sigma[1,1] + SE[1])
        h_m <- Sigma[2,2] / (Sigma[2,2] + SE[2])
        return(c(h_f, h_m))
    } else {
        return(c(NA,NA))
    }    
}

getPlotHeritabilities <- function(cov_dat, fit, type){
    B_dat <- cov_dat$B
    SE_dat <- cov_dat$SE
    N <- cov_dat$N
    Sigma <- getSigma(fit)
    p <- getPi(fit)
    hx <- sapply(1:N, function(i) getH(B_dat[i,], SE_dat[i,], p, Sigma))

    png(sprintf('%s/%s_heritability.png', DATA.DIR, type))
    plot(hx[1,] ~ hx[2,], ylab="h_f", xlab="h_m", 
         main=paste("SNP-level heritabilities (", type, ")", sep=""), col="darkblue")
    lines(x=seq(0,1), y=seq(0,1), add=TRUE)
    dev.off()
    return(hx)
}

getCategories <- function(B, SE, p, Sigma){
    zeros <-c(0,0)
    SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
    p_1 = p[1]*pmnorm(B, zeros, SE_mat)
    p_2 = p[2]*pmnorm(B, zeros, SE_mat + Sigma)
    prob_1 = p_1 / (p_1 + p_2)
    prob_2 = p_2 / (p_1 + p_2)
    category <- ifelse(prob_2 >= 0.8, 2, 1)
    return(category)   
}

overallHeritability <- function(dat, fit){
    cov_dat <- dat$dat
    B_dat <- cov_dat$B
    SE_dat <- cov_dat$SE
    N <- cov_dat$N
    Sigma <- getSigma(fit)
    p <- getPi(fit)
    categories <- sapply(1:N, function(i) getCategories(B_dat[i,], SE_dat[i,], p, Sigma))
    se.p2 <- dat$dat$SE[categories==2,]

    n <- nrow(se.p2)#dat$dat$N
    num_f <- n*(p[2])*Sigma[1,1]
    num_m <- n*(p[2])*Sigma[2,2]

    h_f <- num_f/(num_f + sum(se.p2[,1]))
    h_m <- num_m/(num_m + sum(se.p2[,2]))
    return(list("1"=h_f, "2"=h_m))
}




loadTrait <- function(trait, trait.name){

    print(trait)
    #print(sprintf("Looking at MAF: %s SE: %s", i, j))

    load(sprintf("%s/f_m1_%s_m%s_s%s.RData", DATA.FOLDER, trait, maf, se))
                                        #print(dat$dat$N)
    m1.pi <- getPi(fit1)
    m1.Sigma <- getSigma(fit1)

    print(fit1, pars=c("Sigma", "pi", "Omegacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)                                        #print(m1.pi)
                                        #print(m1.Sigma)
    rg <- getRg(fit1)
    rg.c <- getRgConf(fit1)
    rg.l <- rg.c[[1]]
    rg.u <- rg.c[[2]]


    # assign each SNP to a component, estimate heritability
    h <- overallHeritability(dat, fit1)
    hx <- getPlotHeritabilities(dat$dat, fit1, trait.name)
    next.row <- c(trait, trait.name, dat$dat$N, unlist(m1.pi), unlist(m1.Sigma), rg, rg.l, rg.u, h$'1', h$'2')

    return(next.row)
}

list.traits <- c("48", "20150", "3063", "3064", "4080", "4079", "49", "21001", "50", "whr")
trait.names <- c("WC", "FEV-1", "FVC", "PEF", "BP-s", "BP-d", "HC", "BMI", "height", "whr")
col.labels <- c("trait", "trait.name", "N", "pi[1]", "pi[2]", "Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]", 
    "rg", "rg.l", "rg.u", "h.f", "h.m")

list.out <- lapply(1:length(list.traits), function(i) loadTrait(list.traits[i], trait.names[i]))
df <- do.call(rbind,list.out)
colnames(df) <- col.labels 
write.table(df, file=sprintf("%s/combined_df.txt", DATA.FOLDER))


#se <- 0.4
#list.traits <- c("RH107", "RH160")
#trait.names <- c("Acq_hypo", "asthma")
#write.table(df, file=sprintf("%s/combined_df_bin.txt", DATA.FOLDER))

