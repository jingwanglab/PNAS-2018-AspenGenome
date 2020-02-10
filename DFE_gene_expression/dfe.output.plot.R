#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")

setwd("/proj/b2011141/nobackup/genomic_selection_paper/DFE/plot")
library(gplots)
library(RColorBrewer)

library(RColorBrewer)
colors <- brewer.pal(10,"Paired")[c(1:4,9,10)]
names(colors) <- c("PotraLight", "PotraDark", "PotrsLight", "PotrsDark", "PotriLight", "PotriDark")

#colors <- brewer.pal(9,"Set1")[c(1:9)]

###read.table tremula
#high expression
High_tremula_zero_fold=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremula/zero_fold/summary/tremula.dfe.zero_fold.High.summary",header=T)
High_tremula_utr3=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremula/utr3/summary/tremula.dfe.utr3.High.summary",header=T)
High_tremula_utr5=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremula/utr5/summary/tremula.dfe.utr5.High.summary",header=T)
High_tremula_intron=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremula/intron/summary/tremula.dfe.intron.High.summary",header=T)
#Low expression
Low_tremula_zero_fold=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremula/zero_fold/summary/tremula.dfe.zero_fold.Low.summary",header=T)
Low_tremula_utr3=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremula/utr3/summary/tremula.dfe.utr3.Low.summary",header=T)
Low_tremula_utr5=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremula/utr5/summary/tremula.dfe.utr5.Low.summary",header=T)
Low_tremula_intron=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremula/intron/summary/tremula.dfe.intron.Low.summary",header=T)

###read.table tremuloides
#high expression
High_tremuloides_zero_fold=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremuloides/zero_fold/summary/tremuloides.dfe.zero_fold.High.summary",header=T)
High_tremuloides_utr3=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremuloides/utr3/summary/tremuloides.dfe.utr3.High.summary",header=T)
High_tremuloides_utr5=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremuloides/utr5/summary/tremuloides.dfe.utr5.High.summary",header=T)
High_tremuloides_intron=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/High/tremuloides/intron/summary/tremuloides.dfe.intron.High.summary",header=T)
#Low expression
Low_tremuloides_zero_fold=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremuloides/zero_fold/summary/tremuloides.dfe.zero_fold.Low.summary",header=T)
Low_tremuloides_utr3=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremuloides/utr3/summary/tremuloides.dfe.utr3.Low.summary",header=T)
Low_tremuloides_utr5=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremuloides/utr5/summary/tremuloides.dfe.utr5.Low.summary",header=T)
Low_tremuloides_intron=read.table("/proj/b2011141/nobackup/genomic_selection_paper/DFE/expression/Low/tremuloides/intron/summary/tremuloides.dfe.intron.Low.summary",header=T)


###zero_fold
tremula_zero_fold=as.table(as.matrix(rbind(High_tremula_zero_fold[1,c(9,10,11,12)],Low_tremula_zero_fold[1,c(9,10,11,12)])))
tremula_zero_fold_025=as.table(rbind(apply(High_tremula_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremula_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremula_zero_fold_975=as.table(rbind(apply(High_tremula_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremula_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))


tremuloides_zero_fold=as.table(as.matrix(rbind(High_tremuloides_zero_fold[1,c(9,10,11,12)],Low_tremuloides_zero_fold[1,c(9,10,11,12)])))
tremuloides_zero_fold_025=as.table(rbind(apply(High_tremuloides_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremuloides_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremuloides_zero_fold_975=as.table(rbind(apply(High_tremuloides_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremuloides_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))

###alpha_omega
tremula_zero_fold_alpha=as.table(as.matrix(rbind(High_tremula_zero_fold[1,c(15,16)],Low_tremula_zero_fold[1,c(15,16)])))
tremula_zero_fold_alpha_025=as.table(rbind(apply(High_tremula_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremula_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremula_zero_fold_alpha_975=as.table(rbind(apply(High_tremula_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremula_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.975))))

tremuloides_zero_fold_alpha=as.table(as.matrix(rbind(High_tremuloides_zero_fold[1,c(15,16)],Low_tremuloides_zero_fold[1,c(15,16)])))
tremuloides_zero_fold_alpha_025=as.table(rbind(apply(High_tremuloides_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremuloides_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremuloides_zero_fold_alpha_975=as.table(rbind(apply(High_tremuloides_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremuloides_zero_fold[c(2:201),c(15,16)],2,quantile,probs=c(.975))))


###utr3
tremula_utr3=as.table(as.matrix(rbind(High_tremula_utr3[1,c(9,10,11,12)],Low_tremula_utr3[1,c(9,10,11,12)])))
tremula_utr3_025=as.table(rbind(apply(High_tremula_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremula_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremula_utr3_975=as.table(rbind(apply(High_tremula_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremula_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))


tremuloides_utr3=as.table(as.matrix(rbind(High_tremuloides_utr3[1,c(9,10,11,12)],Low_tremuloides_utr3[1,c(9,10,11,12)])))
tremuloides_utr3_025=as.table(rbind(apply(High_tremuloides_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremuloides_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremuloides_utr3_975=as.table(rbind(apply(High_tremuloides_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremuloides_utr3[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))


###alpha_omega
tremula_utr3_alpha=as.table(as.matrix(rbind(High_tremula_utr3[1,c(15,16)],Low_tremula_utr3[1,c(15,16)])))
tremula_utr3_alpha_025=as.table(rbind(apply(High_tremula_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremula_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremula_utr3_alpha_975=as.table(rbind(apply(High_tremula_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremula_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.975))))

tremuloides_utr3_alpha=as.table(as.matrix(rbind(High_tremuloides_utr3[1,c(15,16)],Low_tremuloides_utr3[1,c(15,16)])))
tremuloides_utr3_alpha_025=as.table(rbind(apply(High_tremuloides_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremuloides_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremuloides_utr3_alpha_975=as.table(rbind(apply(High_tremuloides_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremuloides_utr3[c(2:201),c(15,16)],2,quantile,probs=c(.975))))


###utr5
tremula_utr5=as.table(as.matrix(rbind(High_tremula_utr5[1,c(9,10,11,12)],Low_tremula_utr5[1,c(9,10,11,12)])))
tremula_utr5_025=as.table(rbind(apply(High_tremula_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremula_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremula_utr5_975=as.table(rbind(apply(High_tremula_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremula_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))

tremuloides_utr5=as.table(as.matrix(rbind(High_tremuloides_utr5[1,c(9,10,11,12)],Low_tremuloides_utr5[1,c(9,10,11,12)])))
tremuloides_utr5_025=as.table(rbind(apply(High_tremuloides_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremuloides_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremuloides_utr5_975=as.table(rbind(apply(High_tremuloides_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremuloides_utr5[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))

###alpha_omega
tremula_utr5_alpha=as.table(as.matrix(rbind(High_tremula_utr5[1,c(15,16)],Low_tremula_utr5[1,c(15,16)])))
tremula_utr5_alpha_025=as.table(rbind(apply(High_tremula_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremula_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremula_utr5_alpha_975=as.table(rbind(apply(High_tremula_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremula_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.975))))

tremuloides_utr5_alpha=as.table(as.matrix(rbind(High_tremuloides_utr5[1,c(15,16)],Low_tremuloides_utr5[1,c(15,16)])))
tremuloides_utr5_alpha_025=as.table(rbind(apply(High_tremuloides_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremuloides_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremuloides_utr5_alpha_975=as.table(rbind(apply(High_tremuloides_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremuloides_utr5[c(2:201),c(15,16)],2,quantile,probs=c(.975))))


##Intron
tremula_intron=as.table(as.matrix(rbind(High_tremula_intron[1,c(9,10,11,12)],Low_tremula_intron[1,c(9,10,11,12)])))
tremula_intron_025=as.table(rbind(apply(High_tremula_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremula_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremula_intron_975=as.table(rbind(apply(High_tremula_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremula_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))


tremuloides_intron=as.table(as.matrix(rbind(High_tremuloides_intron[1,c(9,10,11,12)],Low_tremuloides_intron[1,c(9,10,11,12)])))
tremuloides_intron_025=as.table(rbind(apply(High_tremuloides_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(Low_tremuloides_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
tremuloides_intron_975=as.table(rbind(apply(High_tremuloides_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(Low_tremuloides_intron[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))

tremula_intron_alpha=as.table(as.matrix(rbind(High_tremula_intron[1,c(15,16)],Low_tremula_intron[1,c(15,16)])))
tremula_intron_alpha_025=as.table(rbind(apply(High_tremula_intron[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremula_intron[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremula_intron_alpha_975=as.table(rbind(apply(High_tremula_intron[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremula_intron[c(2:201),c(15,16)],2,quantile,probs=c(.975))))

tremuloides_intron_alpha=as.table(as.matrix(rbind(High_tremuloides_intron[1,c(15,16)],Low_tremuloides_intron[1,c(15,16)])))
tremuloides_intron_alpha_025=as.table(rbind(apply(High_tremuloides_intron[c(2:201),c(15,16)],2,quantile,probs=c(.025)),apply(Low_tremuloides_intron[c(2:201),c(15,16)],2,quantile,probs=c(.025))))
tremuloides_intron_alpha_975=as.table(rbind(apply(High_tremuloides_intron[c(2:201),c(15,16)],2,quantile,probs=c(.975)),apply(Low_tremuloides_intron[c(2:201),c(15,16)],2,quantile,probs=c(.975))))


add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }


###functions used to detect significance and plot the significance levels 
sig_plot=function(ci.u,mp){
	#ci.u here is the upper border of the CI, which could be tremula_zero_fold_975
	#mp here is the barplot2 object name, which could be tremula_bp
	y.cord<-rbind(c(ci.u[1,]+0.01),c(apply(ci.u,2,max)+0.02),
          c(apply(ci.u,2,max)+0.02),c(ci.u[2,]+0.01))
	x.cord<-apply(mp,2,function(x) rep(x,each=2))
	sapply(1:4,function(x) lines(x.cord[,x],y.cord[,x]))
}
star_plot=function(ci.u,mp,star){
	x.text<-colMeans(mp)
	y.text<-apply(ci.u,2,max)+0.05
	##wilcox test
	##wilcox.test(High_tremula_zero_fold$Nes_1,Low_tremula_zero_fold$Nes_1)$p.value
	#for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_intron[,i],Low_tremula_intron[,i])$p.value);cat("\n")}
	text(star,x=x.text,y=y.text)
	#text(c("***","***","***","***"),x=x.text,y=y.text)
###*** P<0.0001
###** P<0.01
###* P<0.05

}	


png(filename="tremula_tremuloides.Nes.expression.png",width=6,height=7.5,units='in',res=300)
par(mfrow=c(4,2))
par(mar=c(2,4.5,4,0.5))
###tremula
tremula_bp=barplot2(tremula_zero_fold,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremula_zero_fold_975,ci.l=tremula_zero_fold_025,col=c(colors[2],colors[1]),ylab="Fraction of sites",ylim=c(0,1))
sig_plot(tremula_zero_fold_975,tremula_bp)
star_plot(tremula_zero_fold_975,tremula_bp,c("***","***","***","***"))

mtext("A: 0-fold",side=3,line=0.05,adj=-0.25,font=1.5,cex=1)
tremuloides_bp=barplot2(tremuloides_zero_fold,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremuloides_zero_fold_975,ci.l=tremuloides_zero_fold_025,col=c(colors[4],colors[3]),ylab="Fraction of sites",ylim=c(0,1))
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_zero_fold[,i],Low_tremuloides_zero_fold[,i])$p.value);cat("\n")}
sig_plot(tremuloides_zero_fold_975,tremuloides_bp)
star_plot(tremuloides_zero_fold_975,tremuloides_bp,c("***","***","ns","***"))
##utr3
#tremula
tremula_bp=barplot2(tremula_utr3,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremula_utr3_975,ci.l=tremula_utr3_025,col=c(colors[2],colors[1]),ylab="Fraction of sites",ylim=c(0,1))
mtext("B: 3'UTR",side=3,line=0.05,adj=-0.25,font=1.5,cex=1)
sig_plot(tremula_utr3_975,tremula_bp)
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_utr3[,i],Low_tremula_utr3[,i])$p.value);cat("\n")}
star_plot(tremula_utr3_975,tremula_bp,c("***","***","***","***"))

##tremuloides
tremuloides_bp=barplot2(tremuloides_utr3,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremuloides_utr3_975,ci.l=tremuloides_utr3_025,col=c(colors[4],colors[3]),ylab="Fraction of sites",ylim=c(0,1.05))
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_utr3[,i],Low_tremuloides_utr3[,i])$p.value);cat("\n")}
sig_plot(tremuloides_utr3_975,tremuloides_bp)
star_plot(tremuloides_utr3_975,tremuloides_bp,c("***","***","***","***"))

#utr5
#tremula
tremula_bp=barplot2(tremula_utr5,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremula_utr5_975,ci.l=tremula_utr5_025,col=c(colors[2],colors[1]),ylab="Fraction of sites",ylim=c(0,1))
sig_plot(tremula_utr5_975,tremula_bp)
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_utr5[,i],Low_tremula_utr5[,i])$p.value);cat("\n")}
star_plot(tremula_utr5_975,tremula_bp,c("***","***","***","***"))
mtext("C: 5'UTR",side=3,line=0.05,adj=-0.25,font=1.5,cex=1)
#tremuloides
tremuloides_bp=barplot2(tremuloides_utr5,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremuloides_utr5_975,ci.l=tremuloides_utr5_025,col=c(colors[4],colors[3]),ylab="Fraction of sites",ylim=c(0,1))
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_utr5[,i],Low_tremuloides_utr5[,i])$p.value);cat("\n")}
sig_plot(tremuloides_utr5_975,tremuloides_bp)
star_plot(tremuloides_utr5_975,tremuloides_bp,c("***","***","***","**"))


#tremula_bp=barplot2(tremula_intron,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremula_intron_975,ci.l=tremula_intron_025,col=c(colors[2],colors[1]),ylab="Fraction of sites",xlab=expression(italic(P.tremula)))
#intron
#tremula
tremula_bp=barplot2(tremula_intron,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremula_intron_975,ci.l=tremula_intron_025,col=c(colors[2],colors[1]),ylab="Fraction of sites",ylim=c(0,1))
sig_plot(tremula_intron_975,tremula_bp)
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_intron[,i],Low_tremula_intron[,i])$p.value);cat("\n")}
star_plot(tremula_intron_975,tremula_bp,c("***","***","***","***"))

mtext("D: Intronic",side=3,line=0.05,adj=-0.25,font=1.5,cex=1)
#tremuloides_bp=barplot2(tremuloides_intron,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremuloides_intron_975,ci.l=tremuloides_intron_025,col=c(colors[4],colors[3]),ylab="Fraction of sites",xlab=expression(italic(P.tremuloides)))
#tremuloides
tremuloides_bp=barplot2(tremuloides_intron,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100"))),plot.ci=T,ci.u=tremuloides_intron_975,ci.l=tremuloides_intron_025,col=c(colors[4],colors[3]),ylab="Fraction of sites",ylim=c(0,1))
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_intron[,i],Low_tremuloides_intron[,i])$p.value);cat("\n")}
sig_plot(tremuloides_intron_975,tremuloides_bp)
star_plot(tremuloides_intron_975,tremuloides_bp,c("***","**","ns","***"))

add_legend("top", legend=c(expression(paste("High (",italic(P.tremula),")",sep="")),expression(paste("Low (",italic(P.tremula),")",sep="")),expression(paste("High (",italic(P.tremuloides),")",sep="")),expression(paste("Low (",italic(P.tremuloides),")",sep=""))),bty="n",fill=c(colors[2],colors[1],colors[4],colors[3]),horiz=TRUE)
dev.off()


sig_plot=function(ci.u,mp){
        #ci.u here is the upper border of the CI, which could be tremula_zero_fold_975
        #mp here is the barplot2 object name, which could be tremula_bp
        y.cord<-rbind(c(ci.u[1,]+0.01),c(apply(ci.u,2,max)+0.02),
          c(apply(ci.u,2,max)+0.02),c(ci.u[2,]+0.01))
        x.cord<-apply(mp,2,function(x) rep(x,each=2))
        sapply(1:4,function(x) lines(x.cord[,x],y.cord[,x]))
}


####plot alpha
png(filename="tremula_tremuloides.alpha.expression.png",width=5,height=5,units='in',res=300)
par(mfrow=c(2,1))
par(mar=c(3,4.5,1,1))
#zero_fold_bp=barplot2(as.table(cbind(tremula_zero_fold_alpha[,1],tremuloides_zero_fold_alpha[,1])),beside=T,names.arg=c(expression(italic(P.tremula)),expression(italic(P.tremuloides))),plot.ci=T,ci.u=as.table(cbind(tremula_zero_fold_alpha_975[,1],tremuloides_zero_fold_alpha_975[,1])),ci.l=as.table(cbind(tremula_zero_fold_alpha_025[,1],tremuloides_zero_fold_alpha_025[,1])),col=c("grey40","grey80"),ylab=expression(alpha))
tremula_bp=barplot2(as.table(cbind(tremula_zero_fold_alpha[,1],tremula_utr3_alpha[,1],tremula_utr5_alpha[,1],tremula_intron_alpha[,1])),beside=T,names.arg=c("0-fold","3'UTR","5'UTR","Intronic"),plot.ci=T,ci.u=as.table(cbind(tremula_zero_fold_alpha_975[,1],tremula_utr3_alpha_975[,1],tremula_utr5_alpha_975[,1],tremula_intron_alpha_975[,1])),ci.l=as.table(cbind(tremula_zero_fold_alpha_025[,1],tremula_utr3_alpha_025[,1],tremula_utr5_alpha_025[,1],tremula_intron_alpha_025[,1])),col=c(colors[2],colors[1]),ylab=expression(alpha),ylim=c(0,0.8))
sig_plot(as.table(cbind(tremula_zero_fold_alpha_975[,1],tremula_utr3_alpha_975[,1],tremula_utr5_alpha_975[,1],tremula_intron_alpha_975[,1])),tremula_bp)
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_zero_fold[,i],Low_tremula_zero_fold[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_utr3[,i],Low_tremula_utr3[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_utr5[,i],Low_tremula_utr5[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_intron[,i],Low_tremula_intron[,i])$p.value);cat("\n")}
star_plot(as.table(cbind(tremula_zero_fold_alpha_975[,1],tremula_utr3_alpha_975[,1],tremula_utr5_alpha_975[,1],tremula_intron_alpha_975[,1])),tremula_bp,c("***","***","***","***"))

mtext("(A)",side=3,line=-0.4,adj=-0.2,font=1.5,cex=1)
legend("top", legend=c("High","Low"),bty="n",fill=c(colors[2],colors[1]),horiz=TRUE)
par(mar=c(3,4.5,1,1))
tremuloides_bp=barplot2(as.table(cbind(tremuloides_zero_fold_alpha[,1],tremuloides_utr3_alpha[,1],tremuloides_utr5_alpha[,1],tremuloides_intron_alpha[,1])),beside=T,names.arg=c("0-fold","3'UTR","5'UTR","Intronic"),plot.ci=T,ci.u=as.table(cbind(tremuloides_zero_fold_alpha_975[,1],tremuloides_utr3_alpha_975[,1],tremuloides_utr5_alpha_975[,1],tremuloides_intron_alpha_975[,1])),ci.l=as.table(cbind(tremuloides_zero_fold_alpha_025[,1],tremuloides_utr3_alpha_025[,1],tremuloides_utr5_alpha_025[,1],tremuloides_intron_alpha_025[,1])),col=c(colors[4],colors[3]),ylab=expression(alpha),ylim=c(0,0.8))
sig_plot(as.table(cbind(tremuloides_zero_fold_alpha_975[,1],tremuloides_utr3_alpha_975[,1],tremuloides_utr5_alpha_975[,1],tremuloides_intron_alpha_975[,1])),tremuloides_bp)
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_zero_fold[,i],Low_tremuloides_zero_fold[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_utr3[,i],Low_tremuloides_utr3[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_utr5[,i],Low_tremuloides_utr5[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_intron[,i],Low_tremuloides_intron[,i])$p.value);cat("\n")}
star_plot(as.table(cbind(tremuloides_zero_fold_alpha_975[,1],tremuloides_utr3_alpha_975[,1],tremuloides_utr5_alpha_975[,1],tremuloides_intron_alpha_975[,1])),tremuloides_bp,c("***","***","***","***"))

mtext("(B)",side=3,line=-0.4,adj=-0.2,font=1.5,cex=1)
legend("top", legend=c("High","Low"),bty="n",fill=c(colors[4],colors[3]),horiz=TRUE)
dev.off()



##plot omega
png(filename="tremula_tremuloides.omega.expression.png",width=5,height=5,units='in',res=300)

par(mfrow=c(2,1))
par(mar=c(3,4.5,1,1))
#zero_fold_bp=barplot2(as.table(cbind(tremula_zero_fold_alpha[,1],tremuloides_zero_fold_alpha[,1])),beside=T,names.arg=c(expression(italic(P.tremula)),expression(italic(P.tremuloides))),plot.ci=T,ci.u=as.table(cbind(tremula_zero_fold_alpha_975[,1],tremuloides_zero_fold_alpha_975[,1])),ci.l=as.table(cbind(tremula_zero_fold_alpha_025[,1],tremuloides_zero_fold_alpha_025[,1])),col=c("grey40","grey80"),ylab=expression(alpha))
tremula_bp=barplot2(as.table(cbind(tremula_zero_fold_alpha[,2],tremula_utr3_alpha[,2],tremula_utr5_alpha[,2],tremula_intron_alpha[,2])),beside=T,names.arg=c("0-fold","3'UTR","5'UTR","Intronic"),plot.ci=T,ci.u=as.table(cbind(tremula_zero_fold_alpha_975[,2],tremula_utr3_alpha_975[,2],tremula_utr5_alpha_975[,2],tremula_intron_alpha_975[,2])),ci.l=as.table(cbind(tremula_zero_fold_alpha_025[,2],tremula_utr3_alpha_025[,2],tremula_utr5_alpha_025[,2],tremula_intron_alpha_025[,2])),col=c(colors[2],colors[1]),ylab=expression(omega),ylim=c(0,0.8))
sig_plot(as.table(cbind(tremula_zero_fold_alpha_975[,2],tremula_utr3_alpha_975[,2],tremula_utr5_alpha_975[,2],tremula_intron_alpha_975[,2])),tremula_bp)
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_zero_fold[,i],Low_tremula_zero_fold[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_utr3[,i],Low_tremula_utr3[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_utr5[,i],Low_tremula_utr5[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremula_intron[,i],Low_tremula_intron[,i])$p.value);cat("\n")}
star_plot(as.table(cbind(tremula_zero_fold_alpha_975[,2],tremula_utr3_alpha_975[,2],tremula_utr5_alpha_975[,2],tremula_intron_alpha_975[,2])),tremula_bp,c("***","***","***","***"))
mtext("(A)",side=3,line=-0.4,adj=-0.2,font=1.5,cex=1)
legend("top", legend=c("High","Low"),bty="n",fill=c(colors[2],colors[1]),horiz=TRUE)

par(mar=c(3,4.5,1,1))
tremuloides_bp=barplot2(as.table(cbind(tremuloides_zero_fold_alpha[,2],tremuloides_utr3_alpha[,2],tremuloides_utr5_alpha[,2],tremuloides_intron_alpha[,2])),beside=T,names.arg=c("0-fold","3'UTR","5'UTR","Intronic"),plot.ci=T,ci.u=as.table(cbind(tremuloides_zero_fold_alpha_975[,2],tremuloides_utr3_alpha_975[,2],tremuloides_utr5_alpha_975[,2],tremuloides_intron_alpha_975[,2])),ci.l=as.table(cbind(tremuloides_zero_fold_alpha_025[,2],tremuloides_utr3_alpha_025[,2],tremuloides_utr5_alpha_025[,2],tremuloides_intron_alpha_025[,2])),col=c(colors[4],colors[3]),ylab=expression(omega),ylim=c(0,0.8))

sig_plot(as.table(cbind(tremuloides_zero_fold_alpha_975[,2],tremuloides_utr3_alpha_975[,2],tremuloides_utr5_alpha_975[,2],tremuloides_intron_alpha_975[,2])),tremuloides_bp)
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_zero_fold[,i],Low_tremuloides_zero_fold[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_utr3[,i],Low_tremuloides_utr3[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_utr5[,i],Low_tremuloides_utr5[,i])$p.value);cat("\n")}
for (i in c(9,10,11,12,15,16)){cat(wilcox.test(High_tremuloides_intron[,i],Low_tremuloides_intron[,i])$p.value);cat("\n")}
star_plot(as.table(cbind(tremuloides_zero_fold_alpha_975[,2],tremuloides_utr3_alpha_975[,2],tremuloides_utr5_alpha_975[,2],tremuloides_intron_alpha_975[,2])),tremuloides_bp,c("ns","***","***","***"))

mtext("(B)",side=3,line=-0.4,adj=-0.2,font=1.5,cex=1)
legend("top", legend=c("High","Low"),bty="n",fill=c(colors[4],colors[3]),horiz=TRUE)
dev.off()


