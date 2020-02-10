#! /usr/bin/Rscript --no-save --no-restore

setwd("~/Dropbox/aspen_genome_paper/data/DFE/summary/")
library(gplots)
library(RColorBrewer)

colors <- brewer.pal(10,"Paired")[c(1:10)]

#args=(commandArgs(TRUE))

for (class in c("connectivity","core_gene","eGene","exp_level","exp_variance")) {
  #for (anno in c("0_fold","utr3","utr5","intron","upstream","downstream")) {
  anno="0_fold"
setwd("~/Dropbox/aspen_genome_paper/data/DFE/summary/")
#class=args[1]  ##gene class: connectivity,eGene,core_gene,exp_level,exp_variance
#anno=args[2]  ##annotation,0_fold,utr3,utr5,intron,upstream,downstream

low=read.table(paste("tremula.dfe.",class,".low.",anno,".summary",sep=""),header=T)
high=read.table(paste("tremula.dfe.",class,".high.",anno,".summary",sep=""),header=T)

#####combine dataset
combine=as.table(as.matrix(rbind(low[1,c(9,10,11,12,15,16)],high[1,c(9,10,11,12,15,16)])))
combine_025=as.table(rbind(apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025))))
combine_975=as.table(rbind(apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975))))

###plot
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
  sapply(1:6,function(x) lines(x.cord[,x],y.cord[,x]))
}
star_plot=function(high,low,ci.u,mp){
  x.text<-colMeans(mp)
  y.text<-apply(ci.u,2,max)+0.05
  high=high
  low=low
  ##wilcox test
  #p1=wilcox.test(high$Nes_1,low$Nes_1)$p.value
  p=c()
  for (i in c(9,10,11,12,15,16)){p=c(p,wilcox.test(as.numeric(high[,i]),as.numeric(low[,i]),method="spearman")$p.value)}
                                     #t.test(as.numeric(high[,i]),as.numeric(low[,i]))$p.value)}
                    #wilcox.test(as.numeric(high[,i]),as.numeric(low[,i]),method="spearman")$p.value)}
  p=as.numeric(p)
  star=c()
  for (i in c(1:6)){
    if (p[i] < 0.0001) {star=c(star,"***")}
    else if (p[i] < 0.001 && p[i] > 0.01) {star=c(star,"**")}
    else if (p[i] < 0.01) {star=c(star,"*")}
    else {star=c(star,"")}
  }
  text(star,x=x.text,y=y.text)
  #text(c("***","***","***","***"),x=x.text,y=y.text)
  ###*** P<0.0001
  ###** P<0.001
  ###* P<0.01
}





setwd("~/Dropbox/aspen_genome_paper/data/DFE/plot/")
###plot
pdf(paste("dfe.",class,".",anno,".pdf",sep=""),width=7,height=4)

if (class == "connectivity") {
  col1=colors[1]
  col2=colors[2]
  
par(mar=c(3,5,1,1))
combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2),ylab="Fraction of sites",ylim=c(0,1))
sig_plot(combine_975,combine_bp)
star_plot(high,low,combine_975,combine_bp)

legend("topleft",c("Low-connectivity","High-connectivity"),bty="n",fill=c(col1,col2),horiz=F)
} else if (class == "core_gene") {
  col1=colors[3]
  col2=colors[4]
  
  par(mar=c(3,5,1,1))
  combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2),ylab="Fraction of sites",ylim=c(0,1))
  sig_plot(combine_975,combine_bp)
  star_plot(high,low,combine_975,combine_bp)
  
  legend("top",c("Core","Non-core"),bty="n",fill=c(col1,col2),horiz=F)
} else if (class == "eGene") {
  col1=colors[3]
  col2=colors[4]
  
  par(mar=c(3,5,1,1))
  combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2),ylab="Fraction of sites",ylim=c(0,1))
  sig_plot(combine_975,combine_bp)
  star_plot(low,high,combine_975,combine_bp)
  
  legend("topleft",c("Without eQTLs","With eQTLs"),bty="n",fill=c(col1,col2),horiz=F)
} else if (class == "exp_level") {
  col1=colors[9]
  col2=colors[10]
  
  par(mar=c(3,5,1,1))
  combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2),ylab="Fraction of sites",ylim=c(0,1))
  sig_plot(combine_975,combine_bp)
  star_plot(low,high,combine_975,combine_bp)
  
  legend("topleft",c("Low exp-level","High exp-level"),bty="n",fill=c(col1,col2),horiz=F)
} else if (class == "exp_variance") {
  col1=colors[7]
  col2=colors[8]
  
  par(mar=c(3,5,1,1))
  combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2),ylab="Fraction of sites",ylim=c(0,1))
  sig_plot(combine_975,combine_bp)
  star_plot(low,high,combine_975,combine_bp)
  
  legend("topleft",c("Low exp-variance","High exp-variance"),bty="n",fill=c(col1,col2),horiz=F)
}


dev.off()

  }
}




##################################################################
##output the table
for (class in c("connectivity","eGene","exp_level","core_gene","exp_variance")) {
  for (anno in c("0_fold","intron","utr5","utr3","upstream","downstream")) {
    
    setwd("~/Dropbox/aspen_genome_paper/data/DFE/summary/")
    #class=args[1]  ##gene class: connectivity,eGene,core_gene,exp_level,exp_variance
    #anno=args[2]  ##annotation,0_fold,utr3,utr5,intron,upstream,downstream
    
    low=read.table(paste("tremula.dfe.",class,".low.",anno,".summary",sep=""),header=T)
    high=read.table(paste("tremula.dfe.",class,".high.",anno,".summary",sep=""),header=T)
    
    #####combine dataset
    combine=as.table(as.matrix(rbind(high[1,c(9,10,11,12,15,16)],low[1,c(9,10,11,12,15,16)])))
    combine_025=as.table(rbind(apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025))))
    combine_975=as.table(rbind(apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975))))
    ##output
    cat(class,"\n")
    cat(anno,"\n")
    ##output
    write.table(combine)
    write.table(combine_025)
    write.table(combine_975)
    
    }
}



###Becasue there is quite some different patterns for the efficacy of selection between coding regions and regulatory regions (non-coding), I am wondering
  #it could be caused by different strength of selection acting on the neutral sites (4-fold synonymous sites). Therefore, to be make them in accordance, I 
  #use the same baseline (4-fold synonymous sites in "low" genes from each class) which could likely account for the influence of the choice of neutral sites I guess
######################################################################################
######################################################################################
######################################################################################

setwd("~/Dropbox/aspen_genome_paper/data/DFE/summary_low_neutral/")
library(gplots)
library(RColorBrewer)

colors <- brewer.pal(10,"Paired")[c(1:10)]
colors_blue=brewer.pal(9,"Blues")[c(3,5,7)]
colors_green=brewer.pal(9,"Greens")[c(3,5,7)]
colors_orangle=brewer.pal(9,"Oranges")[c(3,5,7)]
colors_purple=brewer.pal(9,"Purples")[c(3,5,7)]

#args=(commandArgs(TRUE))

for (class in c("connectivity","eGene","exp_level","exp_variance")) {
  #for (anno in c("0_fold","utr3","utr5","intron","upstream","downstream")) {
   anno="0_fold" 
    setwd("~/Dropbox/aspen_genome_paper/data/DFE/summary_low_neutral/")
    #class=args[1]  ##gene class: connectivity,eGene,core_gene,exp_level,exp_variance
    #anno=args[2]  ##annotation,0_fold,utr3,utr5,intron,upstream,downstream
    
    low=read.table(paste("tremula.dfe.",class,".low.",anno,".summary",sep=""),header=T)
    high=read.table(paste("tremula.dfe.",class,".high.",anno,".summary",sep=""),header=T)
    high_neutral=read.table(paste("tremula.dfe.",class,".low_neutral.",anno,".summary",sep=""),header=T)
    
    #####combine dataset
    combine=as.table(as.matrix(rbind(low[1,c(9,10,11,12,15,16)],high[1,c(9,10,11,12,15,16)],high_neutral[1,c(9,10,11,12,15,16)])))
    combine_025=as.table(rbind(apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(high_neutral[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025))))
    combine_975=as.table(rbind(apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(high_neutral[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975))))
    
    ###plot
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
      
      ##here we need two comparisons: (1) low vs. high (2) low vs. high(control)
      ci.u1=ci.u[c(1,2),]
      ci.u2=ci.u[c(1,3),]
      
      mp1=mp[c(1,2),]
      mp2=mp[c(1,3),]
      
      y.cord1<-rbind(c(ci.u1[1,]+0.01),c(apply(ci.u1,2,max)+0.02),
                    c(apply(ci.u1,2,max)+0.02),c(ci.u1[2,]+0.01))
      x.cord1<-apply(mp1,2,function(x) rep(x,each=2))
      sapply(1:6,function(x) lines(x.cord1[,x],y.cord1[,x]))
      
      y.cord2<-rbind(c(ci.u2[1,]+0.08),c(apply(ci.u2,2,max)+0.09),
                     c(apply(ci.u2,2,max)+0.09),c(ci.u2[2,]+0.08))
      x.cord2<-apply(mp2,2,function(x) rep(x,each=2))
      sapply(1:6,function(x) lines(x.cord2[,x],y.cord2[,x]))
      
      }
    
    
    star_plot=function(high,high_neutral,low,ci.u,mp){
      
      ##here we need two comparisons: (1) low vs. high (2) low vs. high(control)
      ci.u1=ci.u[c(1,2),]
      ci.u2=ci.u[c(1,3),]
      
      mp1=mp[c(1,2),]
      mp2=mp[c(1,3),]
      
      
      ### high vs. low
      x.text1<-colMeans(mp1)
      y.text1<-apply(ci.u1,2,max)+0.03
      high=high
      low=low
      ##wilcox test
      #p1=wilcox.test(high$Nes_1,low$Nes_1)$p.value
      p=c()
      ##1. wilcox test
      #for (i in c(9,10,11,12,15,16)){p=c(p,wilcox.test(as.numeric(high[,i]),as.numeric(low[,i]),method="spearman")$p.value)}
      ##2. t.test
      for (i in c(9,10,11,12,15,16)){p=c(p,t.test(as.numeric(high[,i]),as.numeric(low[,i]))$p.value)}
      p=as.numeric(p)
      star=c()
      for (i in c(1:6)){
        if (p[i] < 0.0001) {star=c(star,"***")}
        else if (p[i] < 0.001 && p[i] > 0.01) {star=c(star,"**")}
        else if (p[i] < 0.01) {star=c(star,"*")}
        else {star=c(star,"")}
      }
      text(star,x=x.text1,y=y.text1)
      
      ### high (control) vs. low
      x.text2<-colMeans(mp2)
      y.text2<-apply(ci.u2,2,max)+0.1
      high=high_neutral
      low=low
      ##wilcox test
      #p1=wilcox.test(high$Nes_1,low$Nes_1)$p.value
      p=c()
      #1. wilcox.test
      #for (i in c(9,10,11,12,15,16)){p=c(p,wilcox.test(as.numeric(high[,i]),as.numeric(low[,i]),method="spearman")$p.value)}
      ##2. t.test
      for (i in c(9,10,11,12,15,16)){p=c(p,t.test(as.numeric(high[,i]),as.numeric(low[,i]))$p.value)}
      p=as.numeric(p)
      star=c()
      for (i in c(1:6)){
        if (p[i] < 0.0001) {star=c(star,"***")}
        else if (p[i] < 0.001 && p[i] > 0.01) {star=c(star,"**")}
        else if (p[i] < 0.01) {star=c(star,"*")}
        else {star=c(star,"")}
      }
      text(star,x=x.text2,y=y.text2)
      

      #text(c("***","***","***","***"),x=x.text,y=y.text)
      ###*** P<0.0001
      ###** P<0.001
      ###* P<0.01
    }
    
    
    
    
    
    setwd("~/Dropbox/aspen_genome_paper/data/DFE/plot_low_neutral/")
    ###plot
    pdf(paste("dfe.",class,".",anno,".low_neutral.pdf",sep=""),width=7,height=4)
    
    if (class == "connectivity") {
      col1=colors_blue[1]
      col2=colors_blue[2]
      col3=colors_blue[3]
      
      par(mar=c(3,5,1,1))
      combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2,col3),ylab="Fraction of sites",ylim=c(0,1))
      sig_plot(combine_975,combine_bp)
      star_plot(high,high_neutral,low,combine_975,combine_bp)
      
      legend("topleft",c("Low-connectivity","High-connectivity","High-connectivity (control)"),bty="n",fill=c(col1,col2,col3),horiz=F)
    } else if (class == "eGene") {
    #else if (class == "core_gene") {
    #  col1=colors[3]
    #  col2=colors[4]
    #  
    #  par(mar=c(3,5,1,1))
    #  combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2),ylab="Fraction of sites",ylim=c(0,1))
    #  sig_plot(combine_975,combine_bp)
    #  star_plot(high,low,combine_975,combine_bp)
    #  
    #  legend("top",c("Core","Non-core"),bty="n",fill=c(col1,col2),horiz=TRUE)
    #}
    
      col1=colors_green[1]
      col2=colors_green[2]
      col3=colors_green[3]
      
      
      par(mar=c(3,5,1,1))
      combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2,col3),ylab="Fraction of sites",ylim=c(0,1))
      sig_plot(combine_975,combine_bp)
      star_plot(high,high_neutral,low,combine_975,combine_bp)
      
      legend("topleft",c("Non-eGene","eGene","eGene (control)"),bty="n",fill=c(col1,col2,col3),horiz=F)
    } else if (class == "exp_level") {
      col1=colors_purple[1]
      col2=colors_purple[2]
      col3=colors_purple[3]
      
      par(mar=c(3,5,1,1))
      combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2,col3),ylab="Fraction of sites",ylim=c(0,1))
      sig_plot(combine_975,combine_bp)
      star_plot(high,high_neutral,low,combine_975,combine_bp)
      
      legend("topleft",c("Low exp-level","High exp-level","High exp-level (control)"),bty="n",fill=c(col1,col2,col3),horiz=F)
    } else if (class == "exp_variance") {
      col1=colors_orangle[1]
      col2=colors_orangle[2]
      col3=colors_orangle[3]
      
      par(mar=c(3,5,1,1))
      combine_bp=barplot2(combine,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste("10<",N[e],"s<100")),expression(paste(N[e],"s>100")),expression(alpha),expression(omega)),plot.ci=T,ci.u=combine_975,ci.l=combine_025,col=c(col1,col2,col3),ylab="Fraction of sites",ylim=c(0,1))
      sig_plot(combine_975,combine_bp)
      star_plot(high,high_neutral,low,combine_975,combine_bp)
      
      legend("topleft",c("Low exp-variance","High exp-variance","High exp-variance (control)"),bty="n",fill=c(col1,col2,col3),horiz=F)
    }
    
    
    dev.off()
    
  }
}




##################################################################
##output the table
for (class in c("connectivity","eGene","exp_level","exp_variance")) {
  for (anno in c("0_fold","intron","utr5","utr3","upstream","downstream")) {
    
    setwd("~/Dropbox/aspen_genome_paper/data/DFE/summary_low_neutral/")
    #class=args[1]  ##gene class: connectivity,eGene,core_gene,exp_level,exp_variance
    #anno=args[2]  ##annotation,0_fold,utr3,utr5,intron,upstream,downstream
    
    low=read.table(paste("tremula.dfe.",class,".low.",anno,".summary",sep=""),header=T)
    high=read.table(paste("tremula.dfe.",class,".high.",anno,".summary",sep=""),header=T)
    high_neutral=read.table(paste("tremula.dfe.",class,".low_neutral.",anno,".summary",sep=""),header=T)
    
    #####combine dataset
    combine=as.table(as.matrix(rbind(low[1,c(9,10,11,12,15,16)],high[1,c(9,10,11,12,15,16)],high_neutral[1,c(9,10,11,12,15,16)])))
    combine_025=as.table(rbind(apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025)),apply(high_neutral[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.025))))
    combine_975=as.table(rbind(apply(low[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(high[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975)),apply(high_neutral[c(2:201),c(9,10,11,12,15,16)],2,quantile,probs=c(.975))))
    
    
        ##output
    cat(class,"\n")
    cat(anno,"\n")
    ##output
    write.table(combine)
    write.table(combine_025)
    write.table(combine_975)
    
  }
}




