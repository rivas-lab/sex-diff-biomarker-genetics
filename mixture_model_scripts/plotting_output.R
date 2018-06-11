
require('ggplot2')
require('stringr')
require('tidyr')

summary_files <- list.files(pattern=".txt$")
v2.files <- list.files(pattern="v2.txt$")
v1.files <- setdiff(summary_files, v2.files)
v1.only <- setdiff(v1.files,"ref_data.txt")
s <- lapply(v1.only, function(x) {read.table(x, header=TRUE)})
s.tab <- do.call(rbind, s)
s.tab <- s.tab[s.tab$rg > 0,]

### need to filter out last two - appears they didn't converge ("nucleated ...")

#s.tab$trait <- factor(s.tab$trait, levels=s.tab$trait[order(-s.tab$rg)])
#ggplot(s.tab, aes(x=trait, y=rg))+geom_errorbar(aes(ymin=rg.l, ymax=rg.u)) + geom_point(fill="white")
ref <- read.delim("../multi_run_gwas_input.txt", header=FALSE)
ref2 <- read.delim("../set_traits2.txt", header=FALSE)
ref2$trait.name <- gsub(" ", "_", str_trim(ref2$V3))
ref2$trait <- sapply(ref2$V2, as.character)
ref$trait.name <- gsub(" ", "_", str_trim(paste(ref$V4, ref$V5)))
ref$trait <- sapply(ref$V1, as.character)

ref.total <- rbind(ref[,c("trait", "trait.name")], ref2[,c("trait", "trait.name")])

s.tab$trait <- sapply(s.tab$trait, as.character)

s.tab2 <- merge(ref.total[,c("trait", "trait.name")], s.tab, by="trait")
s.tab2$trait.name <- factor(s.tab2$trait.name, levels=s.tab2$trait.name[order(-s.tab2$rg)])

lit <- read.table("ref_data.txt", header=TRUE)
lit$source <- "literature"
s.tab2$source <- "estimated" 

ref.lit <- merge(lit[,c("rg", "rg.l", "rg.u", "trait.name")], s.tab2[,c("trait.name", "rg", "rg.l", "rg.u")], by="trait.name")

ref.lit2 <- rbind(lit[lit$trait.name %in% ref.lit$trait.name,c("trait.name", "rg", "rg.l", "rg.u", "source")], s.tab2[ s.tab2$trait.name %in% ref.lit$trait.name,c("trait.name", "rg", "rg.l", "rg.u", "source")])

res <- s.tab2[,c("trait.name", "rg", "rg.l", "rg.u", "source")]
ref.lit3 <- rbind(ref.lit2, res)
ref.lit3 <- ref.lit3[!duplicated(ref.lit3),]

# do some fake filtering
to.skip <- sapply(ref.lit3$trait.name, function(x){x <- strsplit(as.character(x), "_", fixed=TRUE)[[1]];  x[[length(x)]] %in% c("volume", 
	"count", "fraction", "concentration", "width", "measure", "mass", "haemoglobin", "crit", "urine", "body", "percentage")})
ref.lit4 <- ref.lit3[!to.skip,]
to.skip2 <- sapply(ref.lit4$trait.name, function(x){x <- strsplit(as.character(x), "_", fixed=TRUE)[[1]];  ifelse(length(x) > 1, x[[(length(x)-1)]] %in% c("leg", "arm", "percentage"), FALSE)})
ref.lit4 <- ref.lit4[!to.skip2,]
ref.lit4$trait.name <- gsub("_", " ", ref.lit4$trait.name)
ref.lit4$trait.name[ref.lit4$trait.name=="FEV-1"] <- "forced expiratory vol"
ref.lit4$trait.name[ref.lit4$trait.name=="PEF"] <- "peak expiratory flow"
ref.lit4$trait.name[ref.lit4$trait.name=="FVC"] <- "forced vital capacity"
ref.lit4$trait.name[ref.lit4$trait.name=="Hand grip strength (right)"] <- "grip strength (r)"

ref.lit4est <- ref.lit4[ref.lit4$source=="estimated",]

ref.lit4$trait.name <- factor(ref.lit4$trait.name, levels=ref.lit4est$trait.name[order(-ref.lit4est$rg)])

#ref.lit2 <- ref.lit %>% gather(source, rg, rg.x:rg.y)
#ref.lit2 <- ref.lit2 %>% gather(source, rg.l, rg.l.x:rg.l.y)
#ref.lit2 <- ref.lit %>% gather(source, rg.u, rg.u.x:rg.u.y)

ggplot(s.tab2, aes(x=trait.name, y=rg))+geom_errorbar(aes(ymin=rg.l, ymax=rg.u)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + geom_point(fill="white")+ylab("Genetic correlation")+xlab("Trait")
ref.lit2$trait.name <- factor(ref.lit2$trait.name, levels=s.tab2$trait.name[order(-s.tab2$rg)])

ggplot(ref.lit4, aes(x=trait.name, y=rg, color=source))+geom_errorbar(aes(ymin=rg.l, ymax=rg.u), width=0.4, size=0.8, position=position_dodge(0.3))+theme(axis.text.x=element_text(angle=90,hjust=1, size=12, face="bold"), legend.title=element_text(size=11, face="bold"), legend.text=element_text(size=11, face="bold"), 
	axis.text.y=element_text(size=11, face="bold"), axis.title.y=element_text(size=12, face="bold")) +ylab("Genetic correlation")+xlab("")+ geom_point(fill="white", position=position_dodge(width=0.1))
# plot the heritabilities
s.tab3 <- s.tab2[,c("trait.name", "h.f", "h.m")]
s.tab3$source <- "estimated"
s.tab3$trait.name <- factor(s.tab3$trait.name, levels=s.tab3$trait.name[order(-s.tab3$h.f)])

h.comb <- rbind(lit[lit$trait.name %in% ref.lit$trait.name,c("trait.name", "h.f", "h.m", "source")], s.tab3[ s.tab2$trait.name %in% ref.lit$trait.name,c("trait.name", "h.f", "h.m", "source")])

s.tab4 <- s.tab3 %>% gather(h.type, h.value, h.f:h.m)
h.comb2 <- h.comb %>% gather(h.type, h.value, h.f:h.m)
h.comb2$h.type.source <- sapply(1:nrow(h.comb2), function(i) paste(h.comb2$h.type[i], h.comb2$source[i], collapse=" "))

h.comb2$trait.name <- factor(h.comb2$trait.name, levels=s.tab3$trait.name[order(-s.tab3$h.f)])
ggplot(h.comb2, aes(x=trait.name, y=h.value, fill=h.type.source))+geom_bar(position=position_dodge(.9), stat="identity")+ylab("Heritability")+xlab("Trait")+scale_fill_manual("legend", values = c("h.f estimated" = "dark blue", "h.f literature" = "light blue", "h.m estimated" = "dark red", "h.m literature" = "red"))

ggplot(h.comb2, aes(x=trait.name, y=h.value, fill=h.type.source))+geom_bar(position=position_dodge(.9), stat="identity")+ylab("Heritability")+xlab("Trait")+  theme(axis.text.x=element_text(angle=90,hjust=1)) +scale_fill_manual("legend", values = c("h.f estimated" = "gray20", "h.f literature" = "gray50", "h.m estimated" = "gray40", "h.m literature" = "gray70"))

#
s.tab5 <- s.tab2[,c("trait.name", "h.f",  "h.auto.f", "h.m", "h.auto.m")]
s.tab5$trait.name <- factor(s.tab5$trait.name, levels=s.tab5$trait.name[order(-s.tab5$h.f)])
colnames(s.tab5) <- c("trait.name", "female - all", "female - autosomal only", "male - all", "male - autosomal only")
s.tab6 <- s.tab5 %>% gather(h.type, h.value, "female - all":"female - autosomal only":"male - all":"male - autosomal only")

ggplot(s.tab6, aes(x=trait.name, y=h.value, fill=h.type))+geom_bar(position=position_dodge(.9), stat="identity")+ylab("Heritability")+xlab("Trait")+  theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_fill_manual("legend", values = c("female - all" = "gray20", "female - autosomal only" = "gray50", "male - all" = "gray40", "male - autosomal only" = "gray70"))
