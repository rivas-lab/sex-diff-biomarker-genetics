#### Model 2 output estimates ####

require('rstan')
require('tidyverse')
require('data.table')
source("src/04_bmm/model_utils.R")
source("src/04_bmm/snp_utils.R")
source("src/04_bmm/heritability_utils.R")
options(stringsAsFactors=FALSE)


# --- FUNCTIONS for grabbing fit output and reformatting --- #
loadDatAnthro <- function(trait){
  load(sprintf("data/vary_priors_%s/m2_fit_2.RData", trait))
  return(fit)
}
loadDatBio <- function(trait, my.dir){
  load(sprintf("%s/m2/f_m2_%s.RData",  my.dir, trait))
  return(fit2)
}

getTraitDf <- function(trait, type="anthro", my.dir="data2"){
    if (type=="anthro"){
        fit <- loadDatAnthro(trait)
    } else {
        fit <- loadDatBio(trait, my.dir)
    }
    summary_df <- data.frame(summary(fit, pars=c("pi", "sigmasq"))$summary)
    summary_df$parameter <- rownames(summary_df)
    summary_df$trait <- trait
    summary_df_long <- summary_df %>% 
        dplyr::select(trait, parameter, "X2.5.", "X50.", "X97.5.", "n_eff", "Rhat") %>%
        rename(ci_l="X2.5.", est="X50.", ci_u="X97.5.") %>%
        as_tibble()
    return(summary_df_long)
}

reformDf <- function(df){
    df %>%
    dplyr::select(trait, parameter, est) %>% 
    tidyr::pivot_wider(id_cols=trait, names_from="parameter", values_from="est") %>%
    rename("pi[0]"="pi[1]",
          "pi[1]"="pi[2]",
          "pi[2]"="pi[3]",
          "pi[3]"="pi[4]",
          "sigmasq[1f]"="sigmasq[1]",
           "sigmasq[2m]"="sigmasq[2]",
           "sigmasq[3f]"="sigmasq[3]",
           "sigmasq[3m]"="sigmasq[4]"
          ) %>%
    mutate_if(is.numeric, ~signif(., digits=4))
}



# --- run on anthro traits --- #
anthro_traits <- c("afr", "tfr", "lfr", "whr")

trait_est <- do.call(rbind, lapply(antrho_traits, getTraitDf))
summary(trait_est$n_eff)
summary(trait_est$Rhat)
trait_reform <- reformDf(trait_est)
trait_reform %>% 
    mutate(trait=case_when(
        trait=="afr" ~ "arm fat ratio",
        trait=="tfr" ~ "trunk fat ratio" ,
        trait=="lfr" ~ "leg fat ratio",
        TRUE ~ trait)) %>%
    write_csv("data/outfiles/m2_anthro_est.csv")


# --- run on bio traits --- #

bio_traits_f <- str_replace_all(list.files(path="data2/m2/", pattern="f*.RData"), "f_m2_|.RData", "")
bio_traits_p <- str_replace_all(list.files(path="data/m2/", pattern="f*.RData"), "f_m2_|.RData", "")
bio_missing_full <- setdiff(bio_traits_p, c(bio_traits_f, anthro_traits))


trait_est_bio_f <- do.call(rbind, lapply(bio_traits_f, getTraitDf, type="bio"))
trait_est_bio_p <- do.call(rbind, lapply(bio_traits_p, getTraitDf, type="bio", my.dir="data"))


# --- clean up, combine, and write out --- #
bio_f <- trait_est_bio_f %>% reformDf()
bio_p <-trait_est_bio_p %>% reformDf()

bio_p %>% 
    anti_join(bio_f, by="trait") %>% 
    filter(!trait %in% anthro_traits) %>%
    bind_rows(bio_f) %>% 
    mutate(trait=str_replace_all(trait, "_", " ")) %>%
    arrange(trait) %>%
    write_csv("data/outfiles/m2_biomarker_est.csv")
