#! /usr/bin/Rscript --no-save --no-restore


setwd("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/plots/LD")

#utr
tremula_utr=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.UTR.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_utr=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.UTR.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_utr=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.UTR.subset.filtering.thin.r2.ld",header=TRUE)
#cds
tremula_cds=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.cds.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_cds=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.cds.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_cds=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.cds.subset.filtering.thin.r2.ld",header=TRUE)
#intron
tremula_intron=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.intron.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_intron=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.intron.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_intron=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.intron.subset.filtering.thin.r2.ld",header=TRUE)
#intergenic
tremula_intergenic=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.intergenic.no_repeats_regulatory.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_intergenic=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.intergenic.no_repeats_regulatory.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_intergenic=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.no_repeats_regulatory.subset.filtering.thin.r2.ld",header=TRUE)
#regulatory
tremula_regulatory=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremula/angsd/plink/subset/tremula.regulatory.subset.filtering.thin.r2.ld",header=TRUE)
tremuloides_regulatory=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/tremuloides/angsd/plink/subset/tremuloides.regulatory.subset.filtering.thin.r2.ld",header=TRUE)
trichocarpa_regulatory=read.table("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/trichocarpa/angsd/plink/subset/trichocarpa.regulatory.subset.filtering.thin.r2.ld",header=TRUE)
#repeats
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
tremula_utr_d=(tremula_utr$BP_B-tremula_utr$BP_A)
tremula_utr_d_order=tremula_utr[order(tremula_utr_d),]
tremula_utr_ld=tremula_utr_d_order$R2

tremula_utr_nlm=nls(tremula_utr_ld~tremula_Er(C_,tremula_utr_d[order(tremula_utr_d)]),start=list(C_=0.01))
tremula_utr_C_<-summary(tremula_utr_nlm)$coefficients[1]
tremula_utr_C_SE<-summary(tremula_utr_nlm)$coefficients[2]


#for P.tremuloides UTR
tremuloides_utr_d=(tremuloides_utr$BP_B-tremuloides_utr$BP_A)
tremuloides_utr_d_order=tremuloides_utr[order(tremuloides_utr_d),]
tremuloides_utr_ld=tremuloides_utr_d_order$R2

tremuloides_utr_nlm=nls(tremuloides_utr_ld~tremuloides_Er(C_,tremuloides_utr_d[order(tremuloides_utr_d)]),start=list(C_=0.01))
tremuloides_utr_C_<-summary(tremuloides_utr_nlm)$coefficients[1]
tremuloides_utr_C_SE<-summary(tremuloides_utr_nlm)$coefficients[2]

#for P.trichocarpa UTR
trichocarpa_utr_d=(trichocarpa_utr$BP_B-trichocarpa_utr$BP_A)
trichocarpa_utr_d_order=trichocarpa_utr[order(trichocarpa_utr_d),]
trichocarpa_utr_ld=trichocarpa_utr_d_order$R2

trichocarpa_utr_nlm=nls(trichocarpa_utr_ld~trichocarpa_Er(C_,trichocarpa_utr_d[order(trichocarpa_utr_d)]),start=list(C_=0.01))
trichocarpa_utr_C_<-summary(trichocarpa_utr_nlm)$coefficients[1]
trichocarpa_utr_C_SE<-summary(trichocarpa_utr_nlm)$coefficients[2]


#cds
#for P.tremula cds
tremula_cds_d=(tremula_cds$BP_B-tremula_cds$BP_A)
tremula_cds_d_order=tremula_cds[order(tremula_cds_d),]
tremula_cds_ld=tremula_cds_d_order$R2

tremula_cds_nlm=nls(tremula_cds_ld~tremula_Er(C_,tremula_cds_d[order(tremula_cds_d)]),start=list(C_=0.01))
tremula_cds_C_<-summary(tremula_cds_nlm)$coefficients[1]
tremula_cds_C_SE<-summary(tremula_cds_nlm)$coefficients[2]

#for P.tremuloides cds
tremuloides_cds_d=(tremuloides_cds$BP_B-tremuloides_cds$BP_A)
tremuloides_cds_d_order=tremuloides_cds[order(tremuloides_cds_d),]
tremuloides_cds_ld=tremuloides_cds_d_order$R2

tremuloides_cds_nlm=nls(tremuloides_cds_ld~tremuloides_Er(C_,tremuloides_cds_d[order(tremuloides_cds_d)]),start=list(C_=0.01))
tremuloides_cds_C_<-summary(tremuloides_cds_nlm)$coefficients[1]
tremuloides_cds_C_SE<-summary(tremuloides_cds_nlm)$coefficients[2]

#for P.trichocarpa cds
trichocarpa_cds_d=(trichocarpa_cds$BP_B-trichocarpa_cds$BP_A)
trichocarpa_cds_d_order=trichocarpa_cds[order(trichocarpa_cds_d),]
trichocarpa_cds_ld=trichocarpa_cds_d_order$R2

trichocarpa_cds_nlm=nls(trichocarpa_cds_ld~trichocarpa_Er(C_,trichocarpa_cds_d[order(trichocarpa_cds_d)]),start=list(C_=0.01))
trichocarpa_cds_C_<-summary(trichocarpa_cds_nlm)$coefficients[1]
trichocarpa_cds_C_SE<-summary(trichocarpa_cds_nlm)$coefficients[2]

#intron

#for P.tremula intron
tremula_intron_d=(tremula_intron$BP_B-tremula_intron$BP_A)
tremula_intron_d_order=tremula_intron[order(tremula_intron_d),]
tremula_intron_ld=tremula_intron_d_order$R2

tremula_intron_nlm=nls(tremula_intron_ld~tremula_Er(C_,tremula_intron_d[order(tremula_intron_d)]),start=list(C_=0.01))
tremula_intron_C_<-summary(tremula_intron_nlm)$coefficients[1]
tremula_intron_C_SE<-summary(tremula_intron_nlm)$coefficients[2]

#for P.tremuloides intron
tremuloides_intron_d=(tremuloides_intron$BP_B-tremuloides_intron$BP_A)
tremuloides_intron_d_order=tremuloides_intron[order(tremuloides_intron_d),]
tremuloides_intron_ld=tremuloides_intron_d_order$R2

tremuloides_intron_nlm=nls(tremuloides_intron_ld~tremuloides_Er(C_,tremuloides_intron_d[order(tremuloides_intron_d)]),start=list(C_=0.01))
tremuloides_intron_C_<-summary(tremuloides_intron_nlm)$coefficients[1]
tremuloides_intron_C_SE<-summary(tremuloides_intron_nlm)$coefficients[2]

#for P.trichocarpa intron
trichocarpa_intron_d=(trichocarpa_intron$BP_B-trichocarpa_intron$BP_A)
trichocarpa_intron_d_order=trichocarpa_intron[order(trichocarpa_intron_d),]
trichocarpa_intron_ld=trichocarpa_intron_d_order$R2

trichocarpa_intron_nlm=nls(trichocarpa_intron_ld~trichocarpa_Er(C_,trichocarpa_intron_d[order(trichocarpa_intron_d)]),start=list(C_=0.01))
trichocarpa_intron_C_<-summary(trichocarpa_intron_nlm)$coefficients[1]
trichocarpa_intron_C_SE<-summary(trichocarpa_intron_nlm)$coefficients[2]

#Intergenic
#for P.tremula Intergenic
tremula_intergenic_d=(tremula_intergenic$BP_B-tremula_intergenic$BP_A)
tremula_intergenic_d_order=tremula_intergenic[order(tremula_intergenic_d),]
tremula_intergenic_ld=tremula_intergenic_d_order$R2

tremula_intergenic_nlm=nls(tremula_intergenic_ld~tremula_Er(C_,tremula_intergenic_d[order(tremula_intergenic_d)]),start=list(C_=0.01))
tremula_intergenic_C_<-summary(tremula_intergenic_nlm)$coefficients[1]
tremula_intergenic_C_SE<-summary(tremula_intergenic_nlm)$coefficients[2]

#for P.tremuloides Intergenic
tremuloides_intergenic_d=(tremuloides_intergenic$BP_B-tremuloides_intergenic$BP_A)
tremuloides_intergenic_d_order=tremuloides_intergenic[order(tremuloides_intergenic_d),]
tremuloides_intergenic_ld=tremuloides_intergenic_d_order$R2

tremuloides_intergenic_nlm=nls(tremuloides_intergenic_ld~tremuloides_Er(C_,tremuloides_intergenic_d[order(tremuloides_intergenic_d)]),start=list(C_=0.01))
tremuloides_intergenic_C_<-summary(tremuloides_intergenic_nlm)$coefficients[1]
tremuloides_intergenic_C_SE<-summary(tremuloides_intergenic_nlm)$coefficients[2]

#for P.trichocarpa Intergenic
trichocarpa_intergenic_d=(trichocarpa_intergenic$BP_B-trichocarpa_intergenic$BP_A)
trichocarpa_intergenic_d_order=trichocarpa_intergenic[order(trichocarpa_intergenic_d),]
trichocarpa_intergenic_ld=trichocarpa_intergenic_d_order$R2

trichocarpa_intergenic_nlm=nls(trichocarpa_intergenic_ld~trichocarpa_Er(C_,trichocarpa_intergenic_d[order(trichocarpa_intergenic_d)]),start=list(C_=0.01))
trichocarpa_intergenic_C_<-summary(trichocarpa_intergenic_nlm)$coefficients[1]
trichocarpa_intergenic_C_SE<-summary(trichocarpa_intergenic_nlm)$coefficients[2]

#Regulatory
#for P.tremula regulatory
tremula_regulatory_d=(tremula_regulatory$BP_B-tremula_regulatory$BP_A)
tremula_regulatory_d_order=tremula_regulatory[order(tremula_regulatory_d),]
tremula_regulatory_ld=tremula_regulatory_d_order$R2

tremula_regulatory_nlm=nls(tremula_regulatory_ld~tremula_Er(C_,tremula_regulatory_d[order(tremula_regulatory_d)]),start=list(C_=0.01))
tremula_regulatory_C_<-summary(tremula_regulatory_nlm)$coefficients[1]
tremula_regulatory_C_SE<-summary(tremula_regulatory_nlm)$coefficients[2]

#for P.tremuloides regulatory
tremuloides_regulatory_d=(tremuloides_regulatory$BP_B-tremuloides_regulatory$BP_A)
tremuloides_regulatory_d_order=tremuloides_regulatory[order(tremuloides_regulatory_d),]
tremuloides_regulatory_ld=tremuloides_regulatory_d_order$R2

tremuloides_regulatory_nlm=nls(tremuloides_regulatory_ld~tremuloides_Er(C_,tremuloides_regulatory_d[order(tremuloides_regulatory_d)]),start=list(C_=0.01))
tremuloides_regulatory_C_<-summary(tremuloides_regulatory_nlm)$coefficients[1]
tremuloides_regulatory_C_SE<-summary(tremuloides_regulatory_nlm)$coefficients[2]

#for P.trichocarpa regulatory
trichocarpa_regulatory_d=(trichocarpa_regulatory$BP_B-trichocarpa_regulatory$BP_A)
trichocarpa_regulatory_d_order=trichocarpa_regulatory[order(trichocarpa_regulatory_d),]
trichocarpa_regulatory_ld=trichocarpa_regulatory_d_order$R2

trichocarpa_regulatory_nlm=nls(trichocarpa_regulatory_ld~trichocarpa_Er(C_,trichocarpa_regulatory_d[order(trichocarpa_regulatory_d)]),start=list(C_=0.01))
trichocarpa_regulatory_C_<-summary(trichocarpa_regulatory_nlm)$coefficients[1]
trichocarpa_regulatory_C_SE<-summary(trichocarpa_regulatory_nlm)$coefficients[2]

#Repeats
tremula_repeats_d=(tremula_repeats$BP_B-tremula_repeats$BP_A)
tremula_repeats_d_order=tremula_repeats[order(tremula_repeats_d),]
tremula_repeats_ld=tremula_repeats_d_order$R2

tremula_repeats_nlm=nls(tremula_repeats_ld~tremula_Er(C_,tremula_repeats_d[order(tremula_repeats_d)]),start=list(C_=0.01))
tremula_repeats_C_<-summary(tremula_repeats_nlm)$coefficients[1]
tremula_repeats_C_SE<-summary(tremula_repeats_nlm)$coefficients[2]

tremuloides_repeats_d=(tremuloides_repeats$BP_B-tremuloides_repeats$BP_A)
tremuloides_repeats_d_order=tremuloides_repeats[order(tremuloides_repeats_d),]
tremuloides_repeats_ld=tremuloides_repeats_d_order$R2

tremuloides_repeats_nlm=nls(tremuloides_repeats_ld~tremuloides_Er(C_,tremuloides_repeats_d[order(tremuloides_repeats_d)]),start=list(C_=0.01))
tremuloides_repeats_C_<-summary(tremuloides_repeats_nlm)$coefficients[1]
tremuloides_repeats_C_SE<-summary(tremuloides_repeats_nlm)$coefficients[2]

trichocarpa_repeats_d=(trichocarpa_repeats$BP_B-trichocarpa_repeats$BP_A)
trichocarpa_repeats_d_order=trichocarpa_repeats[order(trichocarpa_repeats_d),]
trichocarpa_repeats_ld=trichocarpa_repeats_d_order$R2

trichocarpa_repeats_nlm=nls(trichocarpa_repeats_ld~trichocarpa_Er(C_,trichocarpa_repeats_d[order(trichocarpa_repeats_d)]),start=list(C_=0.01))
trichocarpa_repeats_C_<-summary(trichocarpa_repeats_nlm)$coefficients[1]
trichocarpa_repeats_C_SE<-summary(trichocarpa_repeats_nlm)$coefficients[2]


utr_C=cbind(tremula_utr_C_,tremula_utr_C_SE,tremuloides_utr_C_,tremuloides_utr_C_SE,trichocarpa_utr_C_,trichocarpa_utr_C_SE,"UTR")
cds_C=cbind(tremula_cds_C_,tremula_cds_C_SE,tremuloides_cds_C_,tremuloides_cds_C_SE,trichocarpa_cds_C_,trichocarpa_cds_C_SE,"Coding")
intron_C=cbind(tremula_intron_C_,tremula_intron_C_SE,tremuloides_intron_C_,tremuloides_intron_C_SE,trichocarpa_intron_C_,trichocarpa_intron_C_SE,"Intronic")
intergenic_C=cbind(tremula_intergenic_C_,tremula_intergenic_C_SE,tremuloides_intergenic_C_,tremuloides_intergenic_C_SE,trichocarpa_intergenic_C_,trichocarpa_intergenic_C_SE,"Intergenic")
regulatory_C=cbind(tremula_regulatory_C_,tremula_regulatory_C_SE,tremuloides_regulatory_C_,tremuloides_regulatory_C_SE,trichocarpa_regulatory_C_,trichocarpa_regulatory_C_SE,"Regulatory")
repeats_C=cbind(tremula_repeats_C_,tremula_repeats_C_SE,tremuloides_repeats_C_,tremuloides_repeats_C_SE,trichocarpa_repeats_C_,trichocarpa_repeats_C_SE,"Repeats")

total_C=rbind(utr_C,cds_C,intron_C,intergenic_C,regulatory_C,repeats_C)
total_C=as.data.frame(total_C)
names(total_C)=c("tremula_C","tremula_C_SE","tremuloides_C","tremuloides_C_SE","trichocarpa_C","trichocarpa_C_SE","ANNO")

write.table(total_C,file="recombination_rate.se.txt",sep="\t",quote=F,row.names=F,col.names=T)




