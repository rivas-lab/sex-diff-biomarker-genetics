# filterQTdata.R
# E Flynn
# 10/19/2017
# 
# Filter out NAs from QT data

setwd("/scratch/PI/mrivas/users/erflynn/sex_div_gwas")


args = commandArgs(trailingOnly=TRUE)
phe <- args[1]

removeNAs <- function(pheDat){
    # NAs in quantitative data are -9
    print(sprintf("Removed %s missing entries.", table(pheDat$V3 == -9)["TRUE"]))
    return(pheDat[pheDat$V3 != -9,])
}

phe.path <- sprintf("phefiles/%s.phe", phe)

file.copy(phe.path, sprintf("phefiles/old/%s_v0.phe", phe))
pheDat <- read.table(phe.path)
pheDat.f <- removeNAs(pheDat)
write.table(pheDat.f, file=phe.path, col.names=FALSE, row.names=FALSE, quote=FALSE)