#########Ploting DFE figure for different annotation class for P. tremula and P. tremuloides, separately

library(gplots)
library(RColorBrewer)

colors <- brewer.pal(8,"Set2")[c(1:9)]

setwd("~/Dropbox/aspen_genome_paper/data/DFE/tremula_tremuloides/data")

###read.table
tremula_zero_fold=read.table("tremula.dfe.zero_fold.summary",header=T)
tremula_utr3=read.table("tremula.dfe.utr3.summary",header=T)
tremula_utr5=read.table("tremula.dfe.utr5.summary",header=T)
tremula_intron=read.table("tremula.dfe.intron.summary",header=T)
tremula_intergenic=read.table("tremula.dfe.intergenic.summary",header=T)

tremuloides_zero_fold=read.table("tremuloides.dfe.zero_fold.summary",header=T)
tremuloides_utr3=read.table("tremuloides.dfe.utr3.summary",header=T)
tremuloides_utr5=read.table("tremuloides.dfe.utr5.summary",header=T)
tremuloides_intron=read.table("tremuloides.dfe.intron.summary",header=T)
tremuloides_intergenic=read.table("tremuloides.dfe.intergenic.summary",header=T)

###real data
tremula=as.table(as.matrix(rbind(tremula_zero_fold[1,c(9,10,11,12,15,16)],tremula_utr3[1,c(9,10,11,12,15,16)],tremula_utr5[1,c(9,10,11,12,15,16)],tremula_intron[1,c(9,10,11,12,15,16)],tremula_intergenic[1,c(9,10,11,12,15,16)])))
tremula_025=as.table(rbind(apply(tremula_zero_fold[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremula_utr3[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremula_utr5[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremula_intron[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremula_intergenic[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025))))
tremula_975=as.table(rbind(apply(tremula_zero_fold[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremula_utr3[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremula_utr5[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremula_intron[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremula_intergenic[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975))))


tremuloides=as.table(as.matrix(rbind(tremuloides_zero_fold[1,c(9,10,11,12,15,16)],tremuloides_utr3[1,c(9,10,11,12,15,16)],tremuloides_utr5[1,c(9,10,11,12,15,16)],tremuloides_intron[1,c(9,10,11,12,15,16)],tremuloides_intergenic[1,c(9,10,11,12,15,16)])))
tremuloides_025=as.table(rbind(apply(tremuloides_zero_fold[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremuloides_utr3[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremuloides_utr5[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremuloides_intron[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(tremuloides_intergenic[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025))))
tremuloides_975=as.table(rbind(apply(tremuloides_zero_fold[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremuloides_utr3[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremuloides_utr5[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremuloides_intron[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(tremuloides_intergenic[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975))))


pdf("dfe.tremula.pdf",width=7,height=4)
par(mar=c(3,5,1,1))
tremula_bp=barplot2(tremula,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=tremula_975,ci.l=tremula_025,col=c(colors[1],colors[2],colors[3],colors[4],colors[5]),ylim=c(0,1),ylab="Fraction of sites")
legend("top",c("0-fold","3'UTR","5'UTR","Intronic","Intergenic"),bty="n",fill=c(colors[1],colors[2],colors[3],colors[4],colors[5]))
dev.off()

pdf("dfe.tremuloides.pdf",width=7,height=4)
par(mar=c(3,5,1,1))
tremuloides_bp=barplot2(tremuloides,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=tremuloides_975,ci.l=tremuloides_025,col=c(colors[1],colors[2],colors[3],colors[4],colors[5]),ylim=c(0,1),ylab="Fraction of sites")
legend("top",c("0-fold","3'UTR","5'UTR","Intronic","Intergenic"),bty="n",fill=c(colors[1],colors[2],colors[3],colors[4],colors[5]))
dev.off()










