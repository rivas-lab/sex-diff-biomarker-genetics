# computedPhenotypes.R
# E Flynn
# Updated - 11/15/17
#
# Code for computing derived phenotypes

# --- WHR --- #
phe48 <- read.table("phefiles/48.phe") # WC
phe49 <- read.table("phefiles/49.phe") # HC
stopifnot(phe48$V1 == phe49$V1)
whr <- phe48$V3/phe49$V3
pheWhr <- cbind(phe49[,1:2], whr)
write.table(pheWhr, file="phefiles/whr.phe", col.names=FALSE, row.names=FALSE, quote=FALSE)

# --- FEV-1/PVC --- #
phe20150 <- read.table("phefiles/20150.phe") # FEV-1
phe3063 <- read.table("phefiles/3063.phe") # FVC
stopifnot(phe20150$V1 == phe3063$V1)
fp <- phe20150$V3/phe3063$V3
pheFP <- cbind(phe20150[,1:2], fp)
write.table(pheFP, file="phefiles/FEV_FVC.phe", col.names=FALSE, row.names=FALSE, quote=FALSE)

### TODO - FIX!! REMOVE -9 


# --- fat ratios --- #
total <- read.table("/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/ukb11140_20171113_qt/INI23100.phe") # total

leg_l <- read.table("/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/ukb11140_20171113_qt/INI23116.phe") # LEG - L
leg_r <- read.table("/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/ukb11140_20171113_qt/INI23112.phe") # LEG - R
arm_l <- read.table("/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/ukb11140_20171113_qt/INI23124.phe") # ARM - L
arm_r <- read.table("/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/ukb11140_20171113_qt/INI23120.phe") # ARM - R
trunk <- read.table("/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/ukb11140_20171113_qt/INI23128.phe") # TRUNK

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

