# Code for making a Miami plot
# The code is adapted from the Manhattan plot code here:
#   https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

convert_to_bp <- function(dat){
  # this sets up the X axis scale with base pairs
  dat$BPcum <- 0
  s <- 0
  nbp <- c()
  
  for (i in 1:length(unique(dat$CHR))){
    nbp[i] <- max(dat[dat$CHR == i,]$BP)
    dat[dat$CHR == i,"BPcum"] <- dat[dat$CHR == i,"BP"] + s
    s <- s + nbp[i]
  }
  return(dat)
}

prep_miami_dat <-  function(dat1, dat2) {
  dat1.2 <- convert_to_bp(dat1)
  dat2.2 <- convert_to_bp(dat2)
  
  dat1.3 <- dat1.2 %>% mutate(log10P=-log10(P))
  dat2.3 <- dat2.2 %>% mutate(log10P=log10(P))
  
  gwas_dat <- dat1.3 %>% 
    bind_rows(dat2.3) %>% 
    filter(!is.na(P)) %>% 
    mutate(point_grp=((CHR %% 2)==0)) %>%
    mutate(CHR=case_when(
      CHR == 23 ~ "X",  
      CHR == 24  ~"XY",
      CHR == 25 ~ "Y",
      TRUE ~ as.character(CHR))) 
  
  return(gwas_dat)
} 

make_miami_plot <- function(gwas_dat, sig=5*(10**-8)){
  
  # calculate the axis locations
  axis.set <- gwas_dat %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  
  gwas_dat %>% 
    ggplot( aes(x=BPcum, y=log10P)) +
    # alternating colors for variants based on even/odd chromosome
    geom_point(aes(color=as.factor(point_grp)), alpha = 0.3) +
    scale_color_manual(values=c( "gray79", "gray46")) +
    
    # draw lines at y=0 and genome-wide sig for reference
    geom_hline(yintercept = 0, color = "black")+ 
    geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = log10(sig), color = "grey40", linetype = "dashed") +
    # axes
    labs(x = NULL, y = "-log10(p)") + 
    # adjust the X axis scale
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    
    # correct Y axis labels
    scale_y_continuous(breaks=c(-30, -20, -10, 0, 10, 20, 30), 
                       label=c("-30"=">30", "-20"="20", "-10"="10", "0"="0", "10"="10", "20"="20", "30"=">30"))+
    
    # clean up the background and labels
    theme_bw()+
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
   )
}
