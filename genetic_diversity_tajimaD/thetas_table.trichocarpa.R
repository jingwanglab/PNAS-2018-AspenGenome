#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script

library(utils)
library(lattice)
library(RColorBrewer)

args=(commandArgs(TRUE))
anno <- args[1] ##anno here is 5UTR, 3UTR,0_fold,4_fold,intron,regulatory.downstream,regulatory.upstream,intergenic.no_repeats_regulatory
window <- 1000
step <- 1000
trichocarpa=paste("trichocarpa.",anno,".filtering.thetas1kbwindow.1kbsteps.gz.pestPG",sep="")

#set working directory
wd=paste("/proj/b2010014/nobackup/population_genetics/trichocarpa/ANGSD/SFS/folded/anno_type")
setwd(wd)


#set window size and step size
#window <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)
#step <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)

#read table
trichocarpa_thetas <- read.table(trichocarpa)

#set names of columns
names(trichocarpa_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

#normarize theta by numSites
trichocarpa_thetas <- cbind(trichocarpa_thetas, trichocarpa_thetas[,c("tW", "tP", "tF", "tH", "tL")] / trichocarpa_thetas$numSites)
names(trichocarpa_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites", "tW.norm", "tP.norm", "tF.norm", "tH.norm", "tL.norm")


#only selected the sites that have to been 20% left for analysis
minimum_sites=0.1*window
trichocarpa_thetas[which(trichocarpa_thetas$numSites < minimum_sites),c(9,10,11,12,13,15,16,17,18,19)]=rep(NA,10)

trichocarpa_theta=data.frame(cbind(Chr=as.character(trichocarpa_thetas$chr), Pos=trichocarpa_thetas$pos, numSites=trichocarpa_thetas$numSites,tW.norm=trichocarpa_thetas$tW.norm,tP.norm=trichocarpa_thetas$tP.norm,tajD=trichocarpa_thetas$tajD,fulif=trichocarpa_thetas$fulif,fuliD=trichocarpa_thetas$fuliD,fayH=trichocarpa_thetas$fayH,zengsE=trichocarpa_thetas$zengsE))

write.table(trichocarpa_theta, file=paste("trichocarpa_",anno,".",window,"bp",step,"bp",".thetas.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)

