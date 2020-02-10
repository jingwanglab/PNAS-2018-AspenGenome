library(RColorBrewer)


setwd("~/Dropbox/PaperIII/data/gene/v1.1/")

gene_anno=read.csv2("Significant_1615snps.annotation.summary.gene.csv",header=T)
gene=as.character(gene_anno[which(gene_anno$Potra_gene!="Unknown"),"Potra_gene"])
potra_gene=factor(gene)

##1.pie plot of the potra annotation 
#levels(gene_anno$Variant_Classification)=c("IGR","3'Flank","5'Flank","3'UTR","5'UTR","Intron","Splice_Site","Silent","Nonsense_Mutation","Missense_Mutation","Nonstop_Mutation")
potra_anno.table=table(gene_anno$Potra_Variant_Classification)
potra_anno.new=potra_anno.table[c("IGR","3'Flank","5'Flank","3'UTR","5'UTR","Intron","Splice_Site","Silent","Nonsense_Mutation","Missense_Mutation","Nonstop_Mutation")]
cols <- colorRampPalette(brewer.pal(12,"Set3"))(length(potra_anno.table))

percent_potra_anno=round((potra_anno.new/sum(potra_anno.new))*100,digits=2)

labels=names(potra_anno.new)
labels_1=paste(labels,"(",percent_potra_anno,"%)",sep="")

#pdf("snpEFF_maf_pie.pdf")
png(filename="sig_snps.anno_pie.png",width = 5, height = 5, units = 'in', res=300)
par(mfrow=c(1,1))
par(oma=c(0,0,0,0))
par(mgp=c(2,1,0))
par(mar=c(0,3,0,6))

pie(potra_anno.new,col=cols,labels=labels_1,cex=0.75,border=NA,edges=100)
dev.off()

##GO enrichment
source("https://bioconductor.org/biocLite.R")
#biocLite("topGO")
library(topGO)
#library(dplyr)
d=read.table("~/Dropbox/PaperIII/data/gene/topGO/Potra01-protein-in-frame.fa.IPR.GO")
d$V1=as.character(d$V1)
d$V2=as.character(d$V2)
all_go <- lapply(split(d, sub("\\.\\d+$", "", d[, 1])), function(x) unique(x[, 2]))

geneNames=names(all_go)
geneList=factor(as.integer(geneNames %in% potra_gene))
names(geneList)=geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot=annFUN.gene2GO, gene2GO = all_go)
restRes=runTest(GOdata,algorithm="classic",statistic="fisher")
GenTable(GOdata, p.value = restRes, orderBy = "p.value")
gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=20)
write.table(gene_table,file="GO.enrich.txt",sep="\t", quote=F, row.names=F, col.names=T)
