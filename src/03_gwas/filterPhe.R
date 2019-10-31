# filterPhe.R
# E Flynn
# 11/17/2017
#
# Phenotype filtering step for quantitative traits.
# Removes individuals with missing information and/or outside 6sd of mean.

args <- commandArgs(trailingOnly=TRUE)
phe.file <- args[1]
new.phe.file <- args[2]
phe.name <- args[3]
print(sprintf("Phe %s", phe.name))

# remove negatives
phe <- read.table(phe.file)
keep.rows <- (phe[,3] > 0)
phe2 <- phe[keep.rows,]
print(sprintf("Removed %s indiviuals who were missing info.", nrow(phe)-nrow(phe2)))

# remove outside 6 sd
mu <- mean(phe2[,3])
s <- sd(phe2[,3])
cut1 <- mu - 6*s
cut2 <- mu + 6*s
keep.rows2 <- (phe2[,3] > cut1 & phe2[,3] < cut2)
phe3 <- phe2[keep.rows2, ]
print(sprintf("Removed %s indiviuals who were outside 6sd.", nrow(phe2)-nrow(phe3)))

write.table(phe3, file=new.phe.file, row.names=FALSE, col.names=FALSE, quote=FALSE)