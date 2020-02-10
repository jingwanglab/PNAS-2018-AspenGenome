#! /usr/bin/Rscript --no-save --no-restore


setwd("/proj/b2010014/nobackup/private/new_nobackup_area/population_genetics/plots/LD")

total=read.table("recombination_rate.se.txt",header=T)

tremula=as.data.frame(cbind(total$tremula_C,total$tremula_C_SE,as.character(total$ANNO)))
tremula$species="P.tremula"
names(tremula)=c("C","SE","anno","species")

tremuloides=as.data.frame(cbind(total$tremuloides_C,total$tremuloides_C_SE,as.character(total$ANNO)))
tremuloides$species="P.tremuloides"
names(tremuloides)=c("C","SE","anno","species")

trichocarpa=as.data.frame(cbind(total$trichocarpa_C,total$trichocarpa_C_SE,as.character(total$ANNO)))
trichocarpa$species="P.trichocarpa"
names(trichocarpa)=c("C","SE","anno","species")

total_C=rbind(tremula,tremuloides,trichocarpa)
total_C$C=as.numeric(as.character(total_C$C))
total_C$SE=as.numeric(as.character(total_C$SE))

ggplot(total_C,aes(x=anno,y=C,fill=species))+geom_bar(position="dodge",stat="identity")+geom_errorbar(aes(ymin=C-SE,ymax=C+SE),position=position_dodge(0.9),width=.2)+xlab("Annotated features")+ylab(expression(paste("Scaled recombination rate (",rho,")",sep=" ")))+scale_fill_manual(values=c("#1F78B4","#33A02C","#6A3D9A"))

ggsave("recombination_rate.3species.png",width=4,height=4,dpi=300)
