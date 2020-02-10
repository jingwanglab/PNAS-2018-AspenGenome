#! /usr/bin/Rscript --no-save --no-restore

#Here, this R script is used to make LD decay plot, which based on the genotype_r2 produced by vcftools

library(RColorBrewer)
colors <- brewer.pal(10,"Paired")[c(1:4,9,10)]
names(colors) <- c("PotraLight", "PotraDark", "PotrsLight", "PotrsDark", "PotriLight", "PotriDark")

setwd("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/plots/LD")

tremula_utr=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.UTR.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_utr=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.UTR.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_utr=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.UTR.subset.filtering.thin.r2.ld",header=TRUE)

tremula_n=48
tremuloides_n=44
trichocarpa_n=48


#LD Function
tremula_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(tremula_n*(2+C_*d)*(11+C_*d))))
return(res)
}

tremuloides_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(tremuloides_n*(2+C_*d)*(11+C_*d))))
return(res)
}

trichocarpa_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(trichocarpa_n*(2+C_*d)*(11+C_*d))))
return(res)
}


#for P.tremula UTR
tremula_utr_d=(tremula_utr$BP_B-tremula_utr$BP_A)
tremula_utr_d_order=tremula_utr[order(tremula_utr_d),]
tremula_utr_ld=tremula_utr_d_order$R2

tremula_utr_nlm=nls(tremula_utr_ld~tremula_Er(C_,tremula_utr_d[order(tremula_utr_d)]),start=list(C_=0.01))
tremula_utr_C_<-summary(tremula_utr_nlm)$coefficients[1]

#for P.tremuloides UTR
tremuloides_utr_d=(tremuloides_utr$BP_B-tremuloides_utr$BP_A)
tremuloides_utr_d_order=tremuloides_utr[order(tremuloides_utr_d),]
tremuloides_utr_ld=tremuloides_utr_d_order$R2

tremuloides_utr_nlm=nls(tremuloides_utr_ld~tremuloides_Er(C_,tremuloides_utr_d[order(tremuloides_utr_d)]),start=list(C_=0.01))
tremuloides_utr_C_<-summary(tremuloides_utr_nlm)$coefficients[1]

#for P.trichocarpa UTR
trichocarpa_utr_d=(trichocarpa_utr$BP_B-trichocarpa_utr$BP_A)
trichocarpa_utr_d_order=trichocarpa_utr[order(trichocarpa_utr_d),]
trichocarpa_utr_ld=trichocarpa_utr_d_order$R2

trichocarpa_utr_nlm=nls(trichocarpa_utr_ld~trichocarpa_Er(C_,trichocarpa_utr_d[order(trichocarpa_utr_d)]),start=list(C_=0.01))
trichocarpa_utr_C_<-summary(trichocarpa_utr_nlm)$coefficients[1]

#pdf("ld_tremula_tremuloides.pdf")
jpeg("ld_tremula_tremuloides_trichocarpa_utr.jpeg",900,900,pointsize=30)
plot(tremula_utr_d[order(tremula_utr_d)],tremula_utr_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(tremula_utr_d[order(tremula_utr_d)],tremula_Er(tremula_utr_C_,tremula_utr_d[order(tremula_utr_d)]),col=colors[2],lwd=4)
par(new=T)
plot(tremuloides_utr_d[order(tremuloides_utr_d)],tremuloides_utr_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(tremuloides_utr_d[order(tremuloides_utr_d)],tremuloides_Er(tremuloides_utr_C_,tremuloides_utr_d[order(tremuloides_utr_d)]),col=colors[4],lwd=4)
par(new=T)
plot(trichocarpa_utr_d[order(trichocarpa_utr_d)],trichocarpa_utr_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(trichocarpa_utr_d[order(trichocarpa_utr_d)],trichocarpa_Er(trichocarpa_utr_C_,trichocarpa_utr_d[order(trichocarpa_utr_d)]),col=colors[6],lwd=4)

legend(10000,0.8,c("P.tremula","P.tremuloides","P.trichocarpa"),lty=c(1,1,1),lwd=c(4,4,4),col=c(colors[2],colors[4],colors[6]))

dev.off()



