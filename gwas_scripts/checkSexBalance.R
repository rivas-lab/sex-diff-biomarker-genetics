# checkSexBalance.R
# E Flynn
# 11/17/2017
#
# Check the sex balance in the analysis.

BASE.DIR <- '/scratch/PI/mrivas/users/erflynn/sex_div_gwas/'

args <- commandArgs(trailingOnly=TRUE)
phe.file <- args[1]
phe.name <- args[2]

sexLabelPheData <- function(pheDat){
    
    # load phenotype and sex labels
    onesex <- read.table(sprintf('%s/phefiles/onesex.keep', BASE.DIR))
    zerosex <- read.table(sprintf('%s/phefiles/zerosex.keep', BASE.DIR))
    
    pheDat2 <- pheDat[,c(1,3)]
    colnames(pheDat2) <- c("ID", "Phenotype")
    rownames(pheDat2) <- pheDat2$ID
    
    # extract m, f rows
    ids.male <- sapply(unlist(onesex$V1), as.character)
    m <- intersect(ids.male, rownames(pheDat2) )
    men <- pheDat2[m,]
    men$Sex <- 'M'
    ids.female <- sapply(unlist(zerosex$V1), as.character)
    f <- intersect(ids.female, rownames(pheDat2) )
    women <- pheDat2[f,]
    women$Sex <- 'F'
    
    # combine into a full table
    pheFull <- rbind(men, women)
    return(pheFull)
}

# read in, sex label, and count
pheDat <- read.table(phe.file)
pheFull <- sexLabelPheData(pheDat)
sex.counts <- table(pheFull$Sex)
num.f <- sex.counts[["F"]]
num.m <- sex.counts[["M"]]

# output the check
print(sprintf("Phe %s has %s females and %s males", phe.name, num.f, num.m))
if (abs(num.f - num.m) > 1000){
    print(" Please consider sex-balance in this analysis, there is a >1000 subject difference between males and females.")    
}
