#! /usr/bin/Rscript --no-save --no-restore

#Here, this R script is used to make LD decay plot, which based on the genotype_r2 produced by vcftools

library(RColorBrewer)
colors <- brewer.pal(10,"Paired")[c(1:4,9,10)]
names(colors) <- c("PotraLight", "PotraDark", "PotrsLight", "PotrsDark", "PotriLight", "PotriDark")

setwd("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/plots/LD")

tremula_repeats=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.repeats.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_repeats=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.repeats.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_repeats=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.repeats.subset.filtering.thin.r2.ld",header=TRUE)

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
tremula_repeats_d=(tremula_repeats$BP_B-tremula_repeats$BP_A)
tremula_repeats_d_order=tremula_repeats[order(tremula_repeats_d),]
tremula_repeats_ld=tremula_repeats_d_order$R2

tremula_repeats_nlm=nls(tremula_repeats_ld~tremula_Er(C_,tremula_repeats_d[order(tremula_repeats_d)]),start=list(C_=0.01))
tremula_repeats_C_<-summary(tremula_repeats_nlm)$coefficients[1]

#for P.tremuloides UTR
tremuloides_repeats_d=(tremuloides_repeats$BP_B-tremuloides_repeats$BP_A)
tremuloides_repeats_d_order=tremuloides_repeats[order(tremuloides_repeats_d),]
tremuloides_repeats_ld=tremuloides_repeats_d_order$R2

tremuloides_repeats_nlm=nls(tremuloides_repeats_ld~tremuloides_Er(C_,tremuloides_repeats_d[order(tremuloides_repeats_d)]),start=list(C_=0.01))
tremuloides_repeats_C_<-summary(tremuloides_repeats_nlm)$coefficients[1]

#for P.trichocarpa UTR
trichocarpa_repeats_d=(trichocarpa_repeats$BP_B-trichocarpa_repeats$BP_A)
trichocarpa_repeats_d_order=trichocarpa_repeats[order(trichocarpa_repeats_d),]
trichocarpa_repeats_ld=trichocarpa_repeats_d_order$R2

trichocarpa_repeats_nlm=nls(trichocarpa_repeats_ld~trichocarpa_Er(C_,trichocarpa_repeats_d[order(trichocarpa_repeats_d)]),start=list(C_=0.01))
trichocarpa_repeats_C_<-summary(trichocarpa_repeats_nlm)$coefficients[1]

#pdf("ld_tremula_tremuloides.pdf")
jpeg("ld_tremula_tremuloides_trichocarpa_repeats.jpeg",900,900,pointsize=30)
plot(tremula_repeats_d[order(tremula_repeats_d)],tremula_repeats_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(tremula_repeats_d[order(tremula_repeats_d)],tremula_Er(tremula_repeats_C_,tremula_repeats_d[order(tremula_repeats_d)]),col=colors[2],lwd=4)
par(new=T)
plot(tremuloides_repeats_d[order(tremuloides_repeats_d)],tremuloides_repeats_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(tremuloides_repeats_d[order(tremuloides_repeats_d)],tremuloides_Er(tremuloides_repeats_C_,tremuloides_repeats_d[order(tremuloides_repeats_d)]),col=colors[4],lwd=4)
par(new=T)
plot(trichocarpa_repeats_d[order(trichocarpa_repeats_d)],trichocarpa_repeats_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(trichocarpa_repeats_d[order(trichocarpa_repeats_d)],trichocarpa_Er(trichocarpa_repeats_C_,trichocarpa_repeats_d[order(trichocarpa_repeats_d)]),col=colors[6],lwd=4)

legend(10000,0.8,c("P.tremula","P.tremuloides","P.trichocarpa"),lty=c(1,1,1),lwd=c(4,4,4),col=c(colors[2],colors[4],colors[6]))

dev.off()



