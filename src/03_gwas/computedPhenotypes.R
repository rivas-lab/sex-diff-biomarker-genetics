# computedPhenotypes.R
# E Flynn
# Updated - 5/22/2020
#
# Code for computing derived phenotypes

require('tidyverse')

# --- WHR --- #
wc <- read.table("phefiles/INI48.phe")
hc <- read.table("phefiles/INI49.phe")
colnames(wc) <- c("ID", "ID2", "wc")
colnames(hc) <- c("ID", "ID2", "hc")

# remove missing
wc_f <- wc[wc[,3] > 0,]
hc_f <- hc[hc[,3] > 0,]

# remove outside 6 sd
mu <- mean(phe2[,3])
s <- sd(phe2[,3])
cut1 <- mu - 6*s
cut2 <- mu + 6*s
keep.rows2 <- (phe2[,3] > cut1 & phe2[,3] < cut2)
phe3 <- phe2[keep.rows2, ]
print(sprintf("Removed %s indiviuals who were outside 6sd.", nrow(phe2)-nrow(phe3)))


# remove outside 6sd of mean

# calculate


# --- FEV-1/PVC --- #

### TODO - FIX!! REMOVE -9 


# --- fat ratios --- #
total <- read.table("phefiles/INI23100.phe") # total
leg_l <- read.table("phefiles/INI23116.phe")
leg_r <- read.table("phefiles/INI23112.phe")
arm_l <- read.table("phefiles/INI23124.phe")
arm_r <- read.table("phefiles/INI23120.phe")
trunk <- read.table("phefiles/INI23128.phe")

# filter - remove negatives
colnames(leg_l) <- c("ID", "ID2", "leg_l")
colnames(leg_r) <- c("ID", "ID2", "leg_r")
colnames(arm_l) <- c("ID", "ID2", "arm_l")
colnames(arm_r) <- c("ID", "ID2", "arm_r")
colnames(trunk) <- c("ID", "ID2", "trunk")

colnames(total) <- c("ID", "ID2", "total")

leg_l_f <- leg_l[leg_l[,3] > 0,]
leg_r_f <- leg_r[leg_r[,3] > 0,]
arm_l_f <- arm_l[arm_l[,3] > 0,]
arm_r_f <- arm_r[arm_r[,3] > 0,]
trunk_f <- trunk[trunk[,3] > 0,]

total_f <- total[total[,3] >0,]


legs <- merge(leg_l_f[,c(1,3)], leg_r_f[,c(1,3)], by="ID")
arms <- merge(arm_l_f[,c(1,3)], arm_r_f[,c(1,3)], by="ID")
legs.arm <- merge(legs, arms, by="ID")
all <- merge(legs.arm, trunk_f[,c(1,3)], by="ID")

full.df <- merge(all, total_f[,c(1,3)], by="ID")
full.df$leg_sum <- full.df$leg_l + full.df$leg_r
full.df$arm_sum <- full.df$arm_l + full.df$arm_r
full.df$sum <- full.df$leg_sum + full.df$arm_sum + full.df$trunk

 # this is very close to "total" fat mass, but in some cases - it is different
 # max = 4.4, min= -2.05, generally around -0.1 and 0.
 # I am using the sum so it adds to 1 
full.df$leg_p <- full.df$leg_sum/full.df$sum
full.df$arm_p <- full.df$arm_sum/full.df$sum
full.df$trunk_p <- full.df$trunk/full.df$sum

leg_p <- cbind(full.df$ID, full.df$ID, full.df$leg_p)
arm_p <- cbind(full.df$ID, full.df$ID, full.df$arm_p)
trunk_p <- cbind(full.df$ID, full.df$ID, full.df$trunk_p)

write.table(leg_p, file="phefiles/leg_fp.phe", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(arm_p, file="phefiles/arm_fp.phe", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(trunk_p, file="phefiles/trunk_fp.phe", col.names=FALSE, row.names=FALSE, quote=FALSE)

# should we remove outside a certain sd? they did this for phenome paper - outside 5sd
# this would be about 0.1% of individuals

