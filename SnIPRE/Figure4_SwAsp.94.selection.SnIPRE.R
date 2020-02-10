library(data.table)
library(topGO)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ppcor)

setwd("~/Dropbox/aspen_genome_paper/data/SnIPRE_SwAsp94samples/SnIPRE_code_JAGS/SwAsp94_results")

eb=fread("SnIPRE.SwAsp.94samples.eb_results.txt",header=T)
bayes=fread("SnIPRE.SwAsp.94samples.bayesian_results.txt",header=T)
###Using the original expression values in the original files
gene_exp=fread("S2_file_network_eqtl_gene_statistics.tsv",header=T)
gene_exp_all=fread("S2_file_all_gene_statistics.tsv",header=T)

gene_exp_new=gene_exp_all[which(gene_exp_all$gene %in% gene_exp$gene),]
gene_exp$orig_expr_mean=gene_exp_new$expr_mean
gene_exp$orig_expr_var=gene_exp_new$expr_var

##########compare the gene expression characteristics for genes under negative, neutral and positive selection
###selection effect
bayes_pos=bayes[which(bayes$BSnIPRE.class=="pos"),]
bayes_neg=bayes[which(bayes$BSnIPRE.class=="neg"),]
bayes_neut=bayes[which(bayes$BSnIPRE.class=="neut"),]

gene_bayes_pos=gene_exp[which(gene_exp$gene %in% bayes_pos$geneID),]
gene_bayes_neg=gene_exp[which(gene_exp$gene %in% bayes_neg$geneID),]
gene_bayes_neut=gene_exp[which(gene_exp$gene %in% bayes_neut$geneID),]

###constraint effect
bayes_f_pos=bayes[which(bayes$BSnIPRE.f.class=="pos"),]
bayes_f_neg=bayes[which(bayes$BSnIPRE.f.class=="neg"),]
bayes_f_neut=bayes[which(bayes$BSnIPRE.f.class=="neut"),]

gene_bayes_f_pos=gene_exp[which(gene_exp$gene %in% bayes_f_pos$geneID),]
gene_bayes_f_neg=gene_exp[which(gene_exp$gene %in% bayes_f_neg$geneID),]
gene_bayes_f_neut=gene_exp[which(gene_exp$gene %in% bayes_f_neut$geneID),]

boxplot(gene_bayes_f_neg$orig_expr_var,gene_bayes_f_neut$orig_expr_var,outline=F)
boxplot(gene_bayes_f_neg$orig_expr_mean,gene_bayes_f_neut$orig_expr_mean,outline=F)
boxplot(gene_bayes_f_neg$kTotal,gene_bayes_f_neut$kTotal,outline=F)

###########Dividing the genes into negative, neutral and positive selection genes and compare their gene expression characteristics
#selection
bayes_order=bayes[order(bayes$BSnIPRE.est),]
#png("selection_effects.gene.png",width=7,height=5,units='in',res=400)
pdf("selection_effects.gene.pdf",width=7,height=5)
par(mar=c(5,5,1,1))
plot(1:nrow(bayes_order),bayes_order$BSnIPRE.est,type="n",ylim=c(-1,1),xlab="Gene",ylab="Selection effect",cex.axis=1.5,cex.lab=1.5)
##1.negative selected genes
for (n in 1:2109){
segments(x0=n, y0=bayes_order[n,]$BSnIPRE.lbound, x1 = n, y1 = bayes_order[n,]$BSnIPRE.ubound,col="lightblue")
}
for (n in 2110:21911){
segments(x0=n, y0=bayes_order[n,]$BSnIPRE.lbound, x1 = n, y1 = bayes_order[n,]$BSnIPRE.ubound,col="lightgrey")
}
for (n in 21912:22306){
  segments(x0=n, y0=bayes_order[n,]$BSnIPRE.lbound, x1 = n, y1 = bayes_order[n,]$BSnIPRE.ubound,col="red")
}
points(1:nrow(bayes_order),bayes_order$BSnIPRE.est,pch=19,cex=0.3)
abline(h=0,lty=2)
dev.off()

###constraint
bayes_order_f=bayes[order(bayes$BSnIPRE.f),]
#png("constraint_effects.gene.png",width=7,height=5,units='in',res=400)
pdf("constraint_effects.gene.pdf",width=7,height=5)
par(mar=c(5,5,1,1))
plot(1:nrow(bayes_order_f),bayes_order_f$BSnIPRE.f,type="n",ylim=c(0,4),xlab="Gene",ylab="Constraint effect",cex.axis=1.5,cex.lab=1.5)
##1.negative selected genes
for (n in 1:18612){
  segments(x0=n, y0=bayes_order_f[n,]$BSnIPRE.f.lb, x1 = n, y1 = bayes_order_f[n,]$BSnIPRE.f.ub,col="lightblue")
}
for (n in 18613:22302){
  segments(x0=n, y0=bayes_order_f[n,]$BSnIPRE.f.lb, x1 = n, y1 = bayes_order_f[n,]$BSnIPRE.f.ub,col="lightgrey")
}
for (n in 22303:22306){
 # segments(x0=n, y0=bayes_order_f[n,]$BSnIPRE.f.lb, x1 = n, y1 = bayes_order_f[n,]$BSnIPRE.f.ub,col="red")
   segments(x0=n, y0=bayes_order_f[n,]$BSnIPRE.f.lb, x1 = n, y1 = bayes_order_f[n,]$BSnIPRE.f.ub,col="lightgrey")
  
}
points(1:nrow(bayes_order_f),bayes_order_f$BSnIPRE.f,pch=19,cex=0.3)
abline(h=1,lty=2)
dev.off()


####Selection effect
################Making the boxplot 
colors <- brewer.pal(12,"Paired")
###expression variance
#png("selection_effects.gene_var.png",width=3,height=5,units='in',res=400)
pdf("selection_effects.gene_var.pdf",width=3,height=5)
par(mar=c(6,5.5,1,1))
boxplot(gene_bayes_neg$orig_expr_var,gene_bayes_neut$orig_expr_var,gene_bayes_pos$orig_expr_var,outline=F,las=2,
        col=c("lightblue","lightgrey","red"),ylab="Expression variance",names=c("Negative","Neutral","Positive"),
        cex.axis=1.3,cex.lab=1.5,
        ylim=c(0,1.05))
wilcox.test(gene_bayes_neg$orig_expr_var,gene_bayes_neut$orig_expr_var)
wilcox.test(gene_bayes_neut$orig_expr_var,gene_bayes_pos$orig_expr_var)
wilcox.test(gene_bayes_neg$orig_expr_var,gene_bayes_pos$orig_expr_var)

segments(1,0.9,1,0.92)
segments(1,0.92,1.95,0.92)
segments(1.95,0.9,1.95,0.92)
text(1.5,0.94,"***")
segments(2.05,0.9,2.05,0.92)
segments(2.05,0.92,3,0.92)
segments(3,0.9,3,0.92)
text(2.5,0.94,"***")
segments(1,0.96,1,0.98)
segments(1,0.98,3,0.98)
segments(3,0.96,3,0.98)
text(2,1,"***")
dev.off()

###expression mean
#png("selection_effects.gene_mean.png",width=3,height=5,units='in',res=400)
pdf("selection_effects.gene_mean.pdf",width=3,height=5)
par(mar=c(6,5.5,1,1))
boxplot(gene_bayes_neg$orig_expr_mean,gene_bayes_neut$orig_expr_mean,gene_bayes_pos$orig_expr_mean,outline=F,las=2,
        col=c("lightblue","lightgrey","red"),ylab="Expression level",names=c("Negative","Neutral","Positive"),
        cex.axis=1.3,cex.lab=1.5,
        ylim=c(0,11.5))
wilcox.test(gene_bayes_neg$orig_expr_mean,gene_bayes_neut$orig_expr_mean)
wilcox.test(gene_bayes_neut$orig_expr_mean,gene_bayes_pos$orig_expr_mean)
wilcox.test(gene_bayes_neg$orig_expr_mean,gene_bayes_pos$orig_expr_mean)

segments(1,10.2,1,10.4)
segments(1,10.4,1.95,10.4)
segments(1.95,10.2,1.95,10.4)
text(1.5,10.6,"***")
segments(2.05,10.2,2.05,10.4)
segments(2.05,10.4,3,10.4)
segments(3,10.2,3,10.4)
text(2.5,10.6,"***")
segments(1,10.8,1,11)
segments(1,11,3,11)
segments(3,10.8,3,11)
text(2,11.2,"***")

dev.off()
###connectivity
png("selection_effects.connectivity.png",width=3,height=5,units='in',res=400)
par(mar=c(5,5,1,1))
boxplot(gene_bayes_neg$kTotal,gene_bayes_neut$kTotal,gene_bayes_pos$kTotal,outline=F,las=2,
        col=c("lightblue","lightgrey","red"),ylab="Connectivity",names=c("Negative","Neutral","Positive"))
dev.off()

####positive selection are more likely to be the core genes, mean(as.numeric(gene_bayes_pos$is_core_gene))
select_core=read.table(text="Neg Neut Pos core ")
mean(as.numeric(gene_bayes_pos$is_core_gene)) #0.1012658
1-mean(as.numeric(gene_bayes_pos$is_core_gene)) #0.8987342
mean(as.numeric(gene_bayes_neut$is_core_gene)) #0.08160792
1-mean(as.numeric(gene_bayes_neut$is_core_gene)) #0.9183921
mean(as.numeric(gene_bayes_neg$is_core_gene)) #0.06590801
1-mean(as.numeric(gene_bayes_neg$is_core_gene)) #0.934092

#select_core=read.table(text="Negative Neutral Positive 
#                       core 0.0659080 0.08160792 0.1012658
#                       non-core 0.934092 0.9183921 0.8987342",header=T)

#png("selection_effects.core_gene.png",width=3,height=5,units='in',res=400)
pdf("selection_effects.core_gene.pdf",width=3,height=5)
par(mar=c(7,5.5,1,1))
select_core=c(0.0659080,0.08160792,0.1012658)
barplot(select_core,col=c("lightblue","lightgrey","red"),
        #ylab="Prop. of core genes",
        names=c("Negative","Neutral","Postive"),
        cex.axis=1.3,cex.lab=1.5,cex.names=1.3,las=2,
        ylim=c(0,0.135))
title(ylab="Prop. of core genes", line=3.5, cex.lab=1.5)

##*P<0.05
#**P<0.01
#***P<0.001

wilcox.test(as.numeric(gene_bayes_neg$is_core_gene),as.numeric(gene_bayes_neut$is_core_gene))
wilcox.test(as.numeric(gene_bayes_neut$is_core_gene),as.numeric(gene_bayes_pos$is_core_gene))
wilcox.test(as.numeric(gene_bayes_neg$is_core_gene),as.numeric(gene_bayes_pos$is_core_gene))

segments(1,0.105,1,0.11)
segments(1,0.11,1.95,0.11)
segments(1.95,0.105,1.95,0.11)
text(1.5,0.115,"*")
segments(2.05,0.105,2.05,0.11)
segments(2.05,0.11,3,0.11)
segments(3,0.105,3,0.11)
text(2.5,0.115,"ns")
segments(1,0.12,1,0.125)
segments(1,0.125,3,0.125)
segments(3,0.12,3,0.125)
text(2,0.13,"*")
dev.off()

mean(as.numeric(gene_bayes_pos$is_egene)) #0.2708861
mean(as.numeric(gene_bayes_neut$is_egene)) #0.2789617
mean(as.numeric(gene_bayes_neg$is_egene)) #0.2892366


#png("selection_effects.egene.png",width=3,height=5,units='in',res=400)
pdf("selection_effects.egene.pdf",width=3,height=5)
par(mar=c(7,5.5,1,1))
select_egene=c(0.2892366,0.2789617,0.2708861)
barplot(select_egene,col=c("lightblue","lightgrey","red"),
        #ylab="Prop. of genes with eQTLs",
        names=c("Negative","Neutral","Postive"),
        cex.axis=1.3,cex.lab=1.5,cex.names=1.3,las=2,
        ylim=c(0,0.37))
title(ylab="Prop. of genes with eQTLs", line=3.5, cex.lab=1.5)

wilcox.test(as.numeric(gene_bayes_neg$is_egene),as.numeric(gene_bayes_neut$is_egene))
wilcox.test(as.numeric(gene_bayes_neut$is_egene),as.numeric(gene_bayes_pos$is_egene))
wilcox.test(as.numeric(gene_bayes_neg$is_egene),as.numeric(gene_bayes_pos$is_egene))

segments(1,0.3,1,0.31)
segments(1,0.31,1.95,0.31)
segments(1.95,0.3,1.95,0.31)
text(1.5,0.32,"ns")
segments(2.05,0.3,2.05,0.31)
segments(2.05,0.31,3,0.31)
segments(3,0.3,3,0.31)
text(2.5,0.32,"ns")
segments(1,0.33,1,0.34)
segments(1,0.34,3,0.34)
segments(3,0.33,3,0.34)
text(2,0.35,"ns")

dev.off()



###expression variance
#png("constraint_effects.gene_var.png",width=3,height=5,units='in',res=400)
pdf("constraint_effects.gene_var.pdf",width=3,height=5)
par(mar=c(6,5.5,1,1))
boxplot(gene_bayes_f_neg$orig_expr_var,gene_bayes_f_neut$orig_expr_var,outline=F,
        col=c("lightblue","lightgrey"),ylab="Expression variance",
        names=c("Constraint","Neutral"),cex.axis=1.3,cex.lab=1.5,las=2,
        ylim=c(0,0.6))

wilcox.test(as.numeric(gene_bayes_f_neg$orig_expr_var),as.numeric(gene_bayes_f_neut$orig_expr_var))
segments(1,0.57,1,0.58)
segments(1,0.58,1.95,0.58)
segments(1.95,0.57,1.95,0.58)
text(1.5,0.59,"***")

dev.off()

###expression level
#png("constraint_effects.gene_mean.png",width=3,height=5,units='in',res=400)
pdf("constraint_effects.gene_mean.pdf",width=3,height=5)
par(mar=c(6,5.5,1,1))
boxplot(gene_bayes_f_neg$orig_expr_mean,gene_bayes_f_neut$orig_expr_mean,outline=F,
        col=c("lightblue","lightgrey"),ylab="Expression level",
        names=c("Constraint","Neutral"),cex.axis=1.3,cex.lab=1.5,las=2,
        ylim=c(0,11))

wilcox.test(as.numeric(gene_bayes_f_neg$orig_expr_mean),as.numeric(gene_bayes_f_neut$orig_expr_mean))
segments(1,10.4,1,10.6)
segments(1,10.6,1.95,10.6)
segments(1.95,10.4,1.95,10.6)
text(1.5,10.8,"***")

dev.off()

###connectivity
#png("constraint_effects.connectivity.png",width=3,height=5,units='in',res=400)
pdf("constraint_effects.connectivity.pdf",width=3,height=5)
par(mar=c(6,5.5,1,1))
boxplot(gene_bayes_f_neg$kTotal,gene_bayes_f_neut$kTotal,outline=F,
        col=c("lightblue","lightgrey"),ylab="Connectivity",names=c("Constraint","Neutral"),las=2)
dev.off()


mean(as.numeric(gene_bayes_f_neut$is_core_gene)) #0.0498645
mean(as.numeric(gene_bayes_f_neg$is_core_gene))  #0.08650333

wilcox.test(as.numeric(gene_bayes_f_neut$is_core_gene),as.numeric(gene_bayes_f_neg$is_core_gene))

#png("constraint_effects.core_gene.png",width=3,height=5,units='in',res=400)
pdf("constraint_effects.core_gene.pdf",width=3,height=5)
par(mar=c(7,5.5,1,1))
select_core=c(0.08650333,0.0498645)
barplot(select_core,col=c("lightblue","lightgrey"),
        #ylab="Prop. of core genes",
        names=c("Constraint","Neutral"),cex.axis=1.3,cex.lab=1.5,cex.names=1.3,las=2,
        ylim=c(0,0.10))
title(ylab="Prop. of core genes", line=3.5, cex.lab=1.5)
wilcox.test(as.numeric(gene_bayes_f_neg$is_core_gene),as.numeric(gene_bayes_f_neut$is_core_gene))

segments(1,0.09,1,0.093)
segments(1,0.093,1.95,0.093)
segments(1.95,0.09,1.95,0.093)
text(1.5,0.096,"***")
dev.off()


mean(as.numeric(gene_bayes_f_neut$is_egene)) #0.3317073 
mean(as.numeric(gene_bayes_f_neg$is_egene)) #0.2694498

#png("constraint_effects.egene.png",width=3,height=5,units='in',res=400)
pdf("constraint_effects.egene.pdf",width=3,height=5)
par(mar=c(7,5.5,1,1))
select_egene=c(0.2694498,0.3317073)
barplot(select_egene,col=c("lightblue","lightgrey"),
        #ylab="Prop. of genes with eQTLs",
        names=c("Constraint","Neutral"),cex.axis=1.3,cex.lab=1.5,cex.names=1.3,las=2,
        ylim=c(0,0.38))
title(ylab="Prop. of genes with eQTLs", line=3.5, cex.lab=1.5)
wilcox.test(as.numeric(gene_bayes_f_neg$is_egene),as.numeric(gene_bayes_f_neut$is_egene))

segments(1,0.34,1,0.35)
segments(1,0.35,1.95,0.35)
segments(1.95,0.34,1.95,0.35)
text(1.5,0.36,"***")

dev.off()


###Dividing the continuous parameters (kTotal, gene_mean into categories with 10% percentile)
gene_exp=gene_exp %>% mutate(kTotal_quantile=ntile(kTotal,2))
gene_exp=gene_exp %>% mutate(expr_mean_quantile=ntile(expr_mean,2))
gene_exp=gene_exp %>% mutate(expr_var_quantile=ntile(expr_var,2))
gene_exp$BSnIPRE.est=bayes$BSnIPRE.est
gene_exp$BSnIPRE.f=bayes$BSnIPRE.f

###only select the columns that are used to make the plot
gene_exp_sel=select(gene_exp,kTotal_quantile,expr_mean_quantile,expr_var_quantile,is_core_gene,is_egene,BSnIPRE.est,BSnIPRE.f)
###reshape the table to three columns
#1.values of either BSnIPRE.est or BSnIPRE.f
#2.types of values, either BSnIPRE.est or BSnIPRE.f
#2. values of five groups of gene characteristics
#4. types of gene characteristics
gene_exp_sel$is_core_gene=as.numeric(gene_exp_sel$is_core_gene)+1
gene_exp_sel$is_egene=as.numeric(gene_exp_sel$is_egene)+1

gene_exp_sel_melt1=melt(gene_exp_sel, id.vars = c("BSnIPRE.est", "BSnIPRE.f"))
names(gene_exp_sel_melt1)=c("BSnIPRE.est","BSnIPRE.f","gene_charact","gene_group")
gene_exp_sel_melt2=melt(gene_exp_sel_melt1,id.vars=c("gene_charact","gene_group"))

###ggplot
#constrait
constrait=ggplot(gene_exp_sel_melt1, aes(x=factor(gene_group), y=BSnIPRE.f)) + 
  geom_boxplot(aes(fill=factor(gene_group)),outlier.shape=NA)+
  scale_y_continuous(limits=c(0,0.8))
constrait+facet_grid(.~gene_charact)
#selection
selection=ggplot(gene_exp_sel_melt1, aes(x=factor(gene_group), y=BSnIPRE.est)) + 
  geom_boxplot(aes(fill=factor(gene_group)),outlier.shape=NA)+
  scale_y_continuous(limits=c(-0.5,0.5))
selection+facet_grid(.~gene_charact)


#################################################################
###Testing for correlation and partial correlation coefficient between pairs of variables

gene_exp$BSnIPRE.est=bayes$BSnIPRE.est
gene_exp$BSnIPRE.f=bayes$BSnIPRE.f
gene_exp$is_core_gene=as.numeric(gene_exp$is_core_gene)
gene_exp$is_egene=as.numeric(gene_exp$is_egene)

#gene expression level
cor.test(gene_exp$BSnIPRE.est,gene_exp$orig_expr_mean,method="spearman")
pcor.test(gene_exp$BSnIPRE.est,gene_exp$orig_expr_mean,gene_exp[,c("orig_expr_var","kTotal","is_core_gene","is_egene")],method="spearman")

cor.test(gene_exp$BSnIPRE.f,gene_exp$orig_expr_mean,method="spearman")
pcor.test(gene_exp$BSnIPRE.f,gene_exp$orig_expr_mean,gene_exp[,c("orig_expr_var","kTotal","is_core_gene","is_egene")],method="spearman")

#gene expression variance
cor.test(gene_exp$BSnIPRE.est,gene_exp$orig_expr_var,method="spearman")
pcor.test(gene_exp$BSnIPRE.est,gene_exp$orig_expr_var,gene_exp[,c("orig_expr_mean","kTotal","is_core_gene","is_egene")],method="spearman")

cor.test(gene_exp$BSnIPRE.f,gene_exp$orig_expr_var,method="spearman")
pcor.test(gene_exp$BSnIPRE.f,gene_exp$orig_expr_var,gene_exp[,c("orig_expr_mean","kTotal","is_core_gene","is_egene")],method="spearman")

#connectivity
cor.test(gene_exp$BSnIPRE.est,gene_exp$kTotal,method="spearman")
pcor.test(gene_exp$BSnIPRE.est,gene_exp$kTotal,gene_exp[,c("orig_expr_mean","orig_expr_var","is_core_gene","is_egene")],method="spearman")

cor.test(gene_exp$BSnIPRE.f,gene_exp$kTotal,method="spearman")
pcor.test(gene_exp$BSnIPRE.f,gene_exp$kTotal,gene_exp[,c("orig_expr_mean","orig_expr_var","is_core_gene","is_egene")],method="spearman")

#core vs. non-core
cor.test(gene_exp$BSnIPRE.est,gene_exp$is_core_gene,method="spearman")
pcor.test(gene_exp$BSnIPRE.est,gene_exp$is_core_gene,gene_exp[,c("orig_expr_mean","orig_expr_var","kTotal","is_egene")],method="spearman")

cor.test(gene_exp$BSnIPRE.f,gene_exp$is_core_gene,method="spearman")
pcor.test(gene_exp$BSnIPRE.f,gene_exp$is_core_gene,gene_exp[,c("orig_expr_mean","orig_expr_var","kTotal","is_egene")],method="spearman")

#egene vs. non-egene
cor.test(gene_exp$BSnIPRE.est,gene_exp$is_egene,method="spearman")
pcor.test(gene_exp$BSnIPRE.est,gene_exp$is_egene,gene_exp[,c("orig_expr_mean","orig_expr_var","kTotal","is_core_gene")],method="spearman")

cor.test(gene_exp$BSnIPRE.f,gene_exp$is_egene,method="spearman")
pcor.test(gene_exp$BSnIPRE.f,gene_exp$is_egene,gene_exp[,c("orig_expr_mean","orig_expr_var","kTotal","is_core_gene")],method="spearman")






#####GO enrichment for positively selected genes
bayes_pos=bayes[which(bayes$BSnIPRE.class=="pos"),]

d=read.table("../topGO_potra01/Potra01-protein-in-frame.fa.IPR.GO")
d$V1=as.character(d$V1)
d$V2=as.character(d$V2)
all_go <- lapply(split(d, sub("\\.\\d+$", "", d[, 1])), function(x) unique(x[, 2]))

geneNames=names(all_go)
geneList_bayes_pos=factor(as.integer(geneNames %in% bayes_pos$geneID))
names(geneList_bayes_pos)=geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList_bayes_pos,annot=annFUN.gene2GO, gene2GO = all_go)
restRes=runTest(GOdata,algorithm="classic",statistic="fisher")
GenTable(GOdata, p.value = restRes, orderBy = "p.value")
#gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=20)
gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=50)
gene_table_new=gene_table[which(gene_table$Fisher.p<0.05),]
write.table(gene_table_new,file="GO.enrich.SnIPRE.SwAsp.94.pos.txt",sep="\t", quote=F, row.names=F, col.names=T)


#####GO enrichment for genes that are both positively selected and core genes
gene_bayes_pos_core=gene_bayes_pos[which(gene_bayes_pos$is_core_gene=="TRUE"),]

d=read.table("../topGO_potra01/Potra01-protein-in-frame.fa.IPR.GO")
d$V1=as.character(d$V1)
d$V2=as.character(d$V2)
all_go <- lapply(split(d, sub("\\.\\d+$", "", d[, 1])), function(x) unique(x[, 2]))

geneNames=names(all_go)
geneList_bayes_pos_core=factor(as.integer(geneNames %in% gene_bayes_pos_core$gene))
names(geneList_bayes_pos_core)=geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList_bayes_pos_core,annot=annFUN.gene2GO, gene2GO = all_go)
restRes=runTest(GOdata,algorithm="classic",statistic="fisher")
GenTable(GOdata, p.value = restRes, orderBy = "p.value")
gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=50)
gene_table_new=gene_table[which(gene_table$Fisher.p<0.05),]
write.table(gene_table_new,file="GO.enrich.SnIPRE.SwAsp.94.pos.core.txt",sep="\t", quote=F, row.names=F, col.names=T)

cat(as.vector(gene_bayes_pos_core$gene),file="SnIPRE.SwAsp.94.pos.core.gene.txt",sep=",",quote=F)

###
#####GO enrichment for negatively selected genes
bayes_neg=bayes[which(bayes$BSnIPRE.class=="neg"),]

d=read.table("../topGO_potra01/Potra01-protein-in-frame.fa.IPR.GO")
d$V1=as.character(d$V1)
d$V2=as.character(d$V2)
all_go <- lapply(split(d, sub("\\.\\d+$", "", d[, 1])), function(x) unique(x[, 2]))

geneNames=names(all_go)
geneList_bayes_neg=factor(as.integer(geneNames %in% bayes_neg$geneID))
names(geneList_bayes_neg)=geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList_bayes_neg,annot=annFUN.gene2GO, gene2GO = all_go)
restRes=runTest(GOdata,algorithm="classic",statistic="fisher")
GenTable(GOdata, p.value = restRes, orderBy = "p.value")
#gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=20)
gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=100)
gene_table_new=gene_table[which(gene_table$Fisher.p<0.05),]
write.table(gene_table_new,file="GO.enrich.SnIPRE.SwAsp.94.pos.txt",sep="\t", quote=F, row.names=F, col.names=T)





