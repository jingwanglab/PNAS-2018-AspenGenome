#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script

library(utils)
library(lattice)
library(RColorBrewer)

args=(commandArgs(TRUE))
anno <- args[1] ##anno here is 5UTR, 3UTR,0_fold,4_fold,intron,regulatory.downstream,regulatory.upstream,intergenic.no_repeats_regulatory
window <- 1000
step <- 1000
tremuloides=paste("tremuloides.",anno,".filtering.thetas1kbwindow.1kbsteps.gz.pestPG",sep="")

#set working directory
wd=paste("/proj/b2010014/nobackup/population_genetics/tremuloides/ANGSD/SFS/anno_type")
setwd(wd)


#set window size and step size
#window <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)
#step <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)

#read table
tremuloides_thetas <- read.table(tremuloides)

#set names of columns
names(tremuloides_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

#normarize theta by numSites
tremuloides_thetas <- cbind(tremuloides_thetas, tremuloides_thetas[,c("tW", "tP", "tF", "tH", "tL")] / tremuloides_thetas$numSites)
names(tremuloides_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites", "tW.norm", "tP.norm", "tF.norm", "tH.norm", "tL.norm")


#only selected the sites that have to been 20% left for analysis
minimum_sites=0.1*window
tremuloides_thetas[which(tremuloides_thetas$numSites < minimum_sites),c(9,10,11,12,13,15,16,17,18,19)]=rep(NA,10)

tremuloides_theta=data.frame(cbind(Chr=as.character(tremuloides_thetas$chr), Pos=tremuloides_thetas$pos, numSites=tremuloides_thetas$numSites,tW.norm=tremuloides_thetas$tW.norm,tP.norm=tremuloides_thetas$tP.norm,tajD=tremuloides_thetas$tajD,fulif=tremuloides_thetas$fulif,fuliD=tremuloides_thetas$fuliD,fayH=tremuloides_thetas$fayH,zengsE=tremuloides_thetas$zengsE))

write.table(tremuloides_theta, file=paste("tremuloides_",anno,".",window,"bp",step,"bp",".thetas.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)



