require('rstan')
require('tidyverse')
require('loo')

args <- commandArgs(trailingOnly=TRUE)

model <- as.numeric(args[1])
type <- args[2]

if (type=="sim"){
    MY.DIR <- "../../data/vary_priors_sim2/"
} else {
    MY.DIR <- "../../data/vary_priors_8/"
}


getOutput <- function(model, idx){
    print(idx)
    load(sprintf("%s/m%s_fit_%s.RData", MY.DIR, model, idx))
    if (model==1){
       my_summary <- summary(fit, pars=c("Sigma", "pi", "Omegacor", "lp__"))$summary
    } else  {
       my_summary <- summary(fit, pars=c("pi", "sigmasq", "lp__"))$summary
    }
    
    fit_title <- sprintf("m%s_%s", model, idx)

    my_summary2 <- data.frame(my_summary)
    my_summary2$parameter <- rownames(my_summary2)

    my_summary_res <- my_summary2 %>% 
    select(c("parameter", "X2.5.", "X50.", "X97.5.", "n_eff", "Rhat")) %>%
    rename(ci_l="X2.5.", est="X50.", ci_h="X97.5.") %>%
    mutate(fit_params=fit_title)
    return(my_summary_res)
}

num_f <- ifelse(model==1, 13, 8)
output <- lapply(1:num_f, function(i) getOutput(model, i))
my_df <- do.call(rbind, output)
write_csv(my_df, sprintf("../../data/%s_m%s_estimates.csv", type, model))