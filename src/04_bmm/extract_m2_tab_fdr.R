#### Model 2 variant tables + FDR ####

require('rstan')
require('tidyverse')
require('data.table')
source("src/04_bmm/model_utils.R")
source("src/04_bmm/snp_utils.R")
source("src/04_bmm/heritability_utils.R")
options(stringsAsFactors=FALSE)

# ------ Set up the variant annotation table ------ #
variant_tab <- fread('variant_filtering/variant_filter_table.new.tsv.gz', data.table=FALSE)
variant_sm <- variant_tab %>% dplyr::select( ID, POS, REF, ALT, Gene_symbol, maf, Consequence, HGVSp) 


# add in X chromosome info
zs <- fread("data/gwas/ukb24983_v2_hg19.Testosterone_zerosex.genotyped.glm.linear", data.table=FALSE)
Xqc <- fread("data/chr_qc/chrX_qc_table.txt", data.table=FALSE)
XYqc <- fread("data/chr_qc/chrXY_qc_table.txt", data.table=FALSE)

XXY <- rbind(Xqc, XYqc)
XY_g <- zs %>% rename("CHR"="#CHROM") %>% filter(CHR %in% c("X", "XY"))
xDat <- left_join(XY_g, XXY %>% dplyr::select(SNP, A1, A2, MAF), by=c("ID"="SNP"))
xDat2 <- xDat %>% mutate(MAF=1-MAF) %>% dplyr::select(ID, CHR, POS, REF, ALT, MAF) # fix this. 
variant_sm2 <- variant_sm %>% rename(MAF=maf, GENE=Gene_symbol) 


variant_tab <- xDat2 %>% mutate(GENE="", Consequence="", HGVSp="") %>% dplyr::select(colnames(variant_sm2)) %>%
bind_rows(variant_sm2)


loadM2Tab <- function(trait, type="anthro", M2.DIR="data2"){
    if (type=="anthro"){
      df <- read.table(sprintf("../data/vary_priors_%s/sig_snps%s_.txt", trait, trait), sep=" ", header=TRUE); 

    } else {
        df <- read.table(sprintf("../%s/m2/sig_snps%s_.txt", M2.DIR, trait), sep=" ", header=TRUE)
    }

    df %>% 
        rename(ID=SNP, p0=p1, p1=p2, p2=p3, p3=p4) %>%
        mutate(category=category-1) %>%
        mutate(trait=trait)
}

# grab the bio data

m2_tab_bio_f <- do.call(rbind, lapply(bio_traits_f, loadM2Tab, type="bio"))
m2_tab_bio_p <- do.call(rbind, lapply(bio_missing_full, function(x) loadM2Tab(x, type="bio", M2.DIR="data")))
m2_bio <- rbind(m2_tab_bio_f, m2_tab_bio_p)

m2_bio_a <- m2_bio %>% left_join(variant_tab, by=c("ID")) %>%
dplyr::select(trait, ID, CHR, POS,REF,ALT,GENE,MAF,B.f,B.m,SE.f,SE.m,P.f,P.m,p0,p1,p2,p3,category,Consequence,HGVSp)

m2_bio_a2 <- m2_bio_a %>% 
    mutate(POS=as.character(POS)) %>% 
    mutate_if(is.numeric, ~signif(., digits=4)) %>%
    mutate(trait=str_replace_all(trait, "_", " ")) %>%
    arrange(trait, CHR, POS)

m2_bio_a2 %>% write_csv("data/outfiles/m2_bio_snps.csv")


# grab the anthro data
m2_tab_anthro <- do.call(rbind, lapply(anthro_traits, loadM2TabAnthro))

m2_anthro <- m2_tab_anthro %>% left_join(variant_tab, by=c("ID")) %>%
dplyr::select(trait, ID, CHR, POS,REF,ALT,GENE,MAF,B.f,B.m,SE.f,SE.m,P.f,P.m,p0,p1,p2,p3,category,Consequence,HGVSp)


m2_anthro2 <- m2_anthro %>% 
    mutate(POS=as.character(POS)) %>% 
    mutate_if(is.numeric, ~signif(., digits=4)) %>%
    mutate(trait=case_when(
        trait=="afr" ~ "arm fat ratio",
        trait=="tfr" ~ "trunk fat ratio" ,
        trait=="lfr" ~ "leg fat ratio",
        TRUE ~ trait))  %>%    
    arrange(trait, CHR, POS)

m2_anthro2 %>% write_csv("data/outfiles/m2_anthro_snps.csv")


# ---- calculate the FDR ---- #
calcFDRM2 <- function(snps, trait.name, cutoff){

    f_t <- snps %>% filter(trait==trait.name & p1 > cutoff) 
    m_t <- snps %>% filter(trait==trait.name & p2 > cutoff) 
    s_t <- snps %>% filter(trait==trait.name & p3 > cutoff)
    
    f_fdr <- ifelse(length(f_t$p0) > 5, sum(f_t$p0)/(length(f_t$p1)), NA)
    m_fdr <- ifelse(length(m_t$p0) > 5, sum(m_t$p0)/(length(m_t$p2)), NA)
    s_fdr <- ifelse(length(s_t$p0) > 5, sum(s_t$p0)/(length(s_t$p3)), NA)

    return(data.frame(list("trait"=trait.name, "cutoff"=cutoff, "f"=f_fdr, "m"=m_fdr, "s"=s_fdr)))
}


fdr_tabs_anthro <- do.call(rbind, lapply(unique(m2_anthro2$trait), function(trait) 
                   do.call(rbind, lapply(seq(0.5, 0.9, 0.1), 
                                        function(cutoff) calcFDRM2(m2_anthro2, trait, cutoff)))))
fdr_tabs_anthro %>% 
  mutate_if(is.numeric, ~signif(., digits=4)) %>% 
  write_csv("data/outfiles/anthro_fdr.csv")

fdr_tabs_bio <- do.call(rbind, lapply(unique(m2_bio_a2$trait), function(trait) 
                   do.call(rbind, lapply(seq(0.5, 0.9, 0.1), 
                                        function(cutoff) calcFDRM2(m2_bio_a2, trait, cutoff)))))
fdr_tabs_bio %>% 
  mutate_if(is.numeric, ~signif(., digits=4)) %>% 
  write_csv("data/outfiles/bio_fdr.csv")
