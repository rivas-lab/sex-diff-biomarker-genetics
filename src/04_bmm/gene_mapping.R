
require('stringr')
# generate from annotation table - just removes a few columns
annotations <- read.csv("data/reform_ukb_array_ref.csv", colClasses='character')
affy.annote <- annotations[,c("Affy.SNP.ID", "Gene")]
affy.annote.f <- affy.annote[affy.annote$"Affy.SNP.ID" != "---" & affy.annote$"Gene" != " --- ",]
colnames(affy.annote.f)[[1]] <- 'snp'

rs.annote <- annotations[,c("dbSNP.RS.ID", "Gene")]
rs.annote.f <- rs.annote[rs.annote$"dbSNP.RS.ID" != "---" & rs.annote$"Gene" != " --- ",]
colnames(rs.annote.f)[[1]] <- 'snp'

snp.tab <- rbind(affy.annote.f, rs.annote.f)

snp.to.gene <- split(snp.tab$Gene, snp.tab$snp)
snp.to.gene.f <- lapply(snp.to.gene, function(x) unlist(ifelse(length(x) > 1, paste(lapply(x, str_trim), collapse= ","), str_trim(x) )))
snp.gene.table <- data.frame(do.call(rbind, snp.to.gene.f))
colnames(snp.gene.table)[1] <- 'gene' 
snp.gene.table$snp <- names(snp.to.gene.f)
snp.gene.table <- snp.gene.table[order(snp.gene.table$snp),]
rownames(snp.gene.table) <- snp.gene.table$snp
write.table(snp.gene.table, "data/snp_gene_table.txt")