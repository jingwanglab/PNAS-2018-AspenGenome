library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(data.table)
colors <- brewer.pal(9,"Set1")[c(5,2,3)]
colors2 <- brewer.pal(8,"Set2")

setwd("~/Dropbox/aspen_genome_paper/data/ld")

total_table=fread("ld.table.3species.new.txt",header=T)
total_table$anno=factor(total_table$anno,levels=c("Coding","UTR","Intronic","Regulatory","Intergenic","Repeats"))

a=ggplot(total_table,aes(x=Distance,y=r2,group=species,colour=species))+facet_wrap(~anno,ncol=1)+scale_color_manual(values=c("#1F78B4","#33A02C","#6A3D9A"),labels=c(expression(italic("P. tremula")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))+geom_line(aes(colour = species),size=1)+xlab("Distance(bp)")+ylab(expression(r^2))
a
ggsave("ld.3species.anno.a.pdf",height=7,width=5)
b=ggplot(total_table,aes(x=Distance,y=r2,group=anno,colour=anno))+facet_wrap(~species,ncol=1)+scale_color_manual(values=colors2[1:6])+geom_line(aes(colour = anno),size=1)+xlab("Distance(bp)")+ylab(expression(r^2))
b
ggsave("ld.3species.anno.b.pdf",height=7,width=5)






