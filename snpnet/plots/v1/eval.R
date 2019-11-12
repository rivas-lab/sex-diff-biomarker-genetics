fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

###############################################
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))

################################
GBE_ID    <- args[1] # 'HC326'
phe_label <- str_replace(args[2], '_', ' ') # 'MI'
family    <- args[3] # 'binomial'

################################
default_data_dir <- '/oak/stanford/groups/mrivas/projects/degas-risk/final_results/evals/all_z_center_p001_20190805/300'

###############################################
fit_covar_model <- function(in_df, family){
    # fit 
    if(family == 'gaussian'){
        covar_fit <- glm(
            formula = formula('phe ~ 1 + age + sex + PC1 + PC2 + PC3 + PC4'), 
            family = 'gaussian', 
            data = in_df
        )    
    }

    if(family == 'binomial'){ 
        covar_fit <- glm(
            formula = formula('phe ~ 1 + age + sex + PC1 + PC2 + PC3 + PC4'), 
            family = binomial(link="logit"), 
            data = in_df
        )    
    }
    return(covar_fit)
}

compute_OR <- function(df, l_bin, u_bin, middle_df){
    stratified_df <- df %>% 
    filter(l_bin < Percentile, Percentile <= u_bin) %>%
    mutate(PRS = T) %>% select(FID, IID, phe, covar_score, PRS)    

    glmfit <- glm (
        phe ~ covar_score + as.factor(PRS),
        bind_rows(middle_df, stratified_df),
        family=binomial(link="logit")
    )
    
    LOR    <- summary(glmfit)$coefficients[3,1]
    se_LOR <- summary(glmfit)$coefficients[3,2]   
    OR     <- exp(LOR)    
    l_OR   <- exp(LOR - 1.96 * se_LOR)
    u_OR   <- exp(LOR + 1.96 * se_LOR)
        
    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        OR   = OR,
        SE_LOR = se_LOR,
        l_OR = l_OR,
        u_OR = u_OR,
        OR_str = sprintf('%.3f (%.3f-%.3f)', OR, l_OR, u_OR),
#         prevalence = stratified_df %>% select(phe) %>% pull() %>% sum()
    ) %>%
    mutate(OR_str = as.character(OR_str))
}

compute_residuals <- function(df, l_bin, u_bin, middle_df){
    df %>% 
    filter(l_bin < Percentile, Percentile <= u_bin) %>% 
    summarize(
        mean   = mean(phe),
        median = median(phe),
        sd     = sd(phe),
        l_err  = mean - 1.96 * sd,    
        u_err  = mean + 1.96 * sd,
    ) %>%
    mutate(
        l_bin = l_bin,
        u_bin = u_bin,
        mean_str = sprintf('%.3f (%.3f-%.3f)', mean, l_err, u_err)
    ) %>%
    mutate(mean_str = as.character(mean_str))
}

compute_tbl <- function(df, func){
    middle_df <- df %>% 
    filter(0.4 < Percentile, Percentile <= 0.6) %>%
    mutate(PRS = F) %>% 
    select(FID, IID, phe, covar_score, PRS)    

    bind_rows(
        func(df,   0, .01, middle_df),
        func(df, .01, .05, middle_df),
        lapply(2:19, function(x){
            func(df, (x-1)/20, x/20, middle_df)
        }),
        func(df, .95, .99, middle_df),
        func(df, .99, 1, middle_df),
    )    
}

compute_summary_df <- function(df, family){
    if(family == 'gaussian'){
        summary_df <- bind_rows(
            df %>% 
            filter(split == 'test') %>% 
            mutate(FID = IID) %>%
            rename(Percentile = dPRS_p)%>%
            compute_tbl(compute_residuals) %>%
            mutate(PRS_type = 'dPRS'),

            df %>% 
            filter(split == 'test') %>% 
            mutate(FID = IID) %>%
            rename(Percentile = PRS_p)%>%
            compute_tbl(compute_residuals) %>%
            mutate(PRS_type = 'PRS')
        )
    }

    if(family == 'binomial'){ 
        summary_df <- bind_rows(
            df %>% 
            filter(split == 'test') %>% 
            mutate(FID = IID) %>%
            rename(Percentile = dPRS_p)%>%
            compute_tbl(compute_OR) %>%
            mutate(PRS_type = 'dPRS'),

            df %>% 
            filter(split == 'test') %>% 
            mutate(FID = IID) %>%
            rename(Percentile = PRS_p)%>%
            compute_tbl(compute_OR) %>%
            mutate(PRS_type = 'PRS')
        )
    }
    
    return(summary_df)    
}

generate_plot <- function(df, GBE_ID, family, phe_label, middle_line){
    if(family == 'gaussian'){        
        p <- df %>%
        mutate(
            l_err = if_else(sd > 10, NA_real_, l_err),
            u_err = if_else(sd > 10, NA_real_, u_err),
            x_ticks_labels = paste0('[', 100 * l_bin, '% - ', 100 * u_bin, '%]')
        ) %>%
        ggplot(aes(x=reorder(x_ticks_labels, -u_bin), y=mean, color=PRS_type)) +
        geom_point() + 
        geom_errorbar(aes(ymin = l_err, ymax = u_err)) +
        geom_hline(yintercept = middle_line, color='gray')+
        theme_bw() + 
        theme(
            legend.position=c(.15, .9),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)
        ) +
        labs(
            title = sprintf('%s (%s)', GBE_ID, phe_label),
            x = 'The polygenic risk score percentile',
            y = phe_label,
            color = 'Polygenic risk score'
        )
    }

    if(family == 'binomial'){ 
        p <- df %>%
        mutate(
            SE_LOR = if_else(SE_LOR > 10, NA_real_, SE_LOR),
            l_OR = if_else(SE_LOR > 10, NA_real_, l_OR),
            u_OR = if_else(SE_LOR > 10, NA_real_, u_OR),
#             l_SE_OR = exp(log(OR) - SE_LOR),
#             u_SE_OR = exp(log(OR) + SE_LOR),            
            x_ticks_labels = paste0('[', 100 * l_bin, '% - ', 100 * u_bin, '%]')
        ) %>%
        ggplot(aes(x=reorder(x_ticks_labels, -u_bin), y=OR, color=PRS_type)) +
        geom_point() + 
#         geom_errorbar(aes(ymin = l_SE_OR, ymax = u_SE_OR)) + 
        geom_errorbar(aes(ymin = l_OR, ymax = u_OR)) +         
        geom_hline(yintercept = middle_line, color='gray')+
        theme_bw() + 
        theme(
            legend.position=c(.15, .9),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)
        ) +
        labs(
            title = sprintf('%s (%s)', GBE_ID, phe_label),
            x = 'The polygenic risk score percentile',
            y = 'Odds ratio',
            color = 'Polygenic risk score'
        )
    }
        
    return(p)
}

comp_plot <- function(df, GBE_ID, phe_label){
    df %>%
    ggplot(aes(
        x = dPRS, y=PRS,
        color=interaction(dPRS_top,PRS_top,dPRS_bottom,PRS_bottom,sep="-",lex.order=TRUE)
    )) + 
    geom_point(alpha=0.05) + 
    theme_bw() + 
    theme(legend.position = "none") +
    labs(
        title = sprintf('%s (%s)', GBE_ID, phe_label),
        x = 'dPRS', y = 'PRS'
    )

}

comp_summary <- function(df, family, comp_p_thr){
    print(sprintf('The threshold for the top/bottom %d percent', comp_p_thr * 100))    
    df %>% select(PRS, dPRS) %>%
    summarise(
        PRS_top     = quantile(PRS,  probs = comp_p_thr),
        PRS_bottom  = quantile(PRS,  probs = 1 - comp_p_thr),
        dPRS_top    = quantile(dPRS, probs = comp_p_thr),
        dPRS_bottom = quantile(dPRS, probs = 1 - comp_p_thr)
    ) %>% print()
        
    if(family == 'gaussian'){
        print(sprintf('Counts (top)'))
        df %>% count(dPRS_top, PRS_top) %>% print()

        print(sprintf('Counts (bottom)'))
        df %>% count(dPRS_bottom, PRS_bottom) %>% print()
    }
    if(family == 'binomial'){     
        print(sprintf('Counts (top)'))
        df %>% count(dPRS_top, PRS_top, phe) %>% print()

        print(sprintf('Counts (bottom)'))
        df %>% count(dPRS_bottom, PRS_bottom, phe) %>% print()
    }    
}

eval_dPRS_and_PRS <- function(
    GBE_ID, family, phe_label, 
    exts = c('png'),
    comp_p_thr = 0.05, data_dir = default_data_dir
){
    in_df <- fread(file.path(data_dir, paste0(GBE_ID, '.tsv'))) %>% 
    rename(phe = GBE_ID) %>%
    group_by(split) %>%    
    mutate(
        IID = as.integer(IID),
        sex = as.factor(sex),
        phe = phe - if_else(family == 'binomial', 1, 0),
        dPRS_p = rank(-dPRS) / n(),
        PRS_p  = rank(-PRS) / n(),
        dPRS_top = (dPRS_p <= comp_p_thr),
        dPRS_bottom = (1 - comp_p_thr < dPRS_p),
        PRS_top = (PRS_p <= comp_p_thr),
        PRS_bottom = (1 - comp_p_thr < PRS_p)
    ) %>%
    ungroup()

    # fit covariates
    print(sprintf('covariates'))
    covar_fit <- in_df %>%
    filter(split == 'train') %>% 
    select(phe, age, sex, paste0('PC', 1:4)) %>% drop_na() %>%
    fit_covar_model(family)
    print(summary(covar_fit))

    # compute covariate score and the residuals
    test_df <- in_df %>% 
    mutate(
        covar_score = predict(covar_fit, newdata = in_df),
        residual = as.numeric(phe) - covar_score
    ) %>%
    filter(split == 'test')
    
    print('we have test_df')
    test_df %>% colnames() %>% print()
     
    # compare dPRS vs. PRS    
    test_df %>% comp_summary(family, comp_p_thr)
    
    comp <- test_df %>% comp_plot(GBE_ID, phe_label)

    for(ext in exts){
        comp + ggsave(file.path(data_dir, paste0(GBE_ID, '.comp.', ext)))
    }

    # summarize the dPRS/PRS performance
    print(colnames(test_df))
    
    summary_df <- test_df %>% compute_summary_df(family)
    if(family == 'gaussian'){
        middle_line <- (test_df %>% select(phe) %>% pull() %>% mean())        
    }
    if(family == 'binomial'){     
        middle_line <- 1
    }        
    
    p <- summary_df %>% generate_plot(GBE_ID, family, phe_label, middle_line)
    
    summary_df %>% fwrite(file.path(data_dir, paste0(GBE_ID, '.eval.tsv')), sep='\t')

    for(ext in exts){
        p + ggsave(file.path(data_dir, paste0(GBE_ID, '.eval.', ext)))
    }
    
    return(summary_df)
}

###############################################

df <- eval_dPRS_and_PRS(
    GBE_ID = GBE_ID,
    family = family,
    phe_label = phe_label
)
