#! /usr/bin/Rscript --no-save --no-restore

#Here, this R script is used to make LD decay plot, which based on the genotype_r2 produced by vcftools

library(RColorBrewer)
colors <- brewer.pal(10,"Paired")[c(1:4,9,10)]
names(colors) <- c("PotraLight", "PotraDark", "PotrsLight", "PotrsDark", "PotriLight", "PotriDark")

setwd("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/plots/LD")

tremula_intergenic=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.intergenic.no_repeats_regulatory.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_intergenic=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.intergenic.no_repeats_regulatory.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_intergenic=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.no_repeats_regulatory.subset.filtering.thin.r2.ld",header=TRUE)

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


#for P.tremula Intergenic
tremula_intergenic_d=(tremula_intergenic$BP_B-tremula_intergenic$BP_A)
tremula_intergenic_d_order=tremula_intergenic[order(tremula_intergenic_d),]
tremula_intergenic_ld=tremula_intergenic_d_order$R2

tremula_intergenic_nlm=nls(tremula_intergenic_ld~tremula_Er(C_,tremula_intergenic_d[order(tremula_intergenic_d)]),start=list(C_=0.01))
tremula_intergenic_C_<-summary(tremula_intergenic_nlm)$coefficients[1]

#for P.tremuloides Intergenic
tremuloides_intergenic_d=(tremuloides_intergenic$BP_B-tremuloides_intergenic$BP_A)
tremuloides_intergenic_d_order=tremuloides_intergenic[order(tremuloides_intergenic_d),]
tremuloides_intergenic_ld=tremuloides_intergenic_d_order$R2

tremuloides_intergenic_nlm=nls(tremuloides_intergenic_ld~tremuloides_Er(C_,tremuloides_intergenic_d[order(tremuloides_intergenic_d)]),start=list(C_=0.01))
tremuloides_intergenic_C_<-summary(tremuloides_intergenic_nlm)$coefficients[1]

#for P.trichocarpa Intergenic
trichocarpa_intergenic_d=(trichocarpa_intergenic$BP_B-trichocarpa_intergenic$BP_A)
trichocarpa_intergenic_d_order=trichocarpa_intergenic[order(trichocarpa_intergenic_d),]
trichocarpa_intergenic_ld=trichocarpa_intergenic_d_order$R2

trichocarpa_intergenic_nlm=nls(trichocarpa_intergenic_ld~trichocarpa_Er(C_,trichocarpa_intergenic_d[order(trichocarpa_intergenic_d)]),start=list(C_=0.01))
trichocarpa_intergenic_C_<-summary(trichocarpa_intergenic_nlm)$coefficients[1]

#pdf("ld_tremula_tremuloides.pdf")
jpeg("ld_tremula_tremuloides_trichocarpa_intergenic.jpeg",900,900,pointsize=30)
plot(tremula_intergenic_d[order(tremula_intergenic_d)],tremula_intergenic_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(tremula_intergenic_d[order(tremula_intergenic_d)],tremula_Er(tremula_intergenic_C_,tremula_intergenic_d[order(tremula_intergenic_d)]),col=colors[2],lwd=4)
par(new=T)
plot(tremuloides_intergenic_d[order(tremuloides_intergenic_d)],tremuloides_intergenic_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(tremuloides_intergenic_d[order(tremuloides_intergenic_d)],tremuloides_Er(tremuloides_intergenic_C_,tremuloides_intergenic_d[order(tremuloides_intergenic_d)]),col=colors[4],lwd=4)
par(new=T)
plot(trichocarpa_intergenic_d[order(trichocarpa_intergenic_d)],trichocarpa_intergenic_ld,type="n",cex=.5,pch=19,col="grey",xlab="Distance (bp)",ylab=expression(r^2))
lines(trichocarpa_intergenic_d[order(trichocarpa_intergenic_d)],trichocarpa_Er(trichocarpa_intergenic_C_,trichocarpa_intergenic_d[order(trichocarpa_intergenic_d)]),col=colors[6],lwd=4)

legend(10000,0.8,c("P.tremula","P.tremuloides","P.trichocarpa"),lty=c(1,1,1),lwd=c(4,4,4),col=c(colors[2],colors[4],colors[6]))

dev.off()



