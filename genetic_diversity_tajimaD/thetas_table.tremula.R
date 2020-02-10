#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script

library(utils)
library(lattice)
library(RColorBrewer)

args=(commandArgs(TRUE))
anno <- args[1] ##anno here is 5UTR, 3UTR,0_fold,4_fold,intron,regulatory.downstream,regulatory.upstream,intergenic.no_repeats_regulatory
window <- 1000
step <- 1000
tremula=paste("tremula.",anno,".filtering.thetas1kbwindow.1kbsteps.gz.pestPG",sep="")

#set working directory
wd=paste("/proj/b2010014/nobackup/population_genetics/tremula/ANGSD/SFS/anno_type")
setwd(wd)

args <- commandArgs(TRUE)
if (length(args) != 1) {
    message("Usage: plotThetas.R tremula tremuloides!")
    quit("yes")
}

#set window size and step size
#window <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)
#step <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)

#read table
tremula_thetas <- read.table(tremula)

#set names of columns
names(tremula_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

#normarize theta by numSites
tremula_thetas <- cbind(tremula_thetas, tremula_thetas[,c("tW", "tP", "tF", "tH", "tL")] / tremula_thetas$numSites)
names(tremula_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites", "tW.norm", "tP.norm", "tF.norm", "tH.norm", "tL.norm")


#only selected the sites that have to been 20% left for analysis
minimum_sites=0.1*window
tremula_thetas[which(tremula_thetas$numSites < minimum_sites),c(9,10,11,12,13,15,16,17,18,19)]=rep(NA,10)

tremula_theta=data.frame(cbind(Chr=as.character(tremula_thetas$chr), Pos=tremula_thetas$pos, numSites=tremula_thetas$numSites,tW.norm=tremula_thetas$tW.norm,tP.norm=tremula_thetas$tP.norm,tajD=tremula_thetas$tajD,fulif=tremula_thetas$fulif,fuliD=tremula_thetas$fuliD,fayH=tremula_thetas$fayH,zengsE=tremula_thetas$zengsE))

write.table(tremula_theta, file=paste("tremula_",anno,".",window,"bp",step,"bp",".thetas.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)



