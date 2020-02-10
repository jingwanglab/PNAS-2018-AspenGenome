#! /usr/bin/Rscript --no-save --no-restore

.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.3")
#.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.5")

##########################################################################################################
#####Maim aim, output file as the input file format needed for SnIPRE
##columns are:geneID, PR (polymorphic 0_fold), FR (fixed 0_fold), PS (polymorphic 4_fold), FS (fixed 4_fold), Tsil (total 4_fold sites), Trepl (total 0_fold_sites), nout (1, number of outgroup), npop (number of samples, 94)


library(data.table)
library(dplyr)
library(reshape2)


args=(commandArgs(TRUE))
wd=args[1]
set(wd)
gene=fread("SwAsp_94samples.SnIPRE.input.txt",header=T)

T_4_fold=sum(gene$total_4_fold) ##total number of 4_fold synonymous sites
T_0_fold=sum(gene$total_0_fold) ##total number of 0_fold nonsynonymous sites
nout=1
npop=94
###geneID,PR,FR,PS,FS,Tsil,Trepl,nout,npop
out=gene[,c(1,2,3,5,6),with=F]
names(out)=c("geneID","PR","FR","PS","FS")
out$Tsil=rep(T_4_fold,nrow(out))
out$Trepl=rep(T_0_fold,nrow(out))
out$nout=rep(nout,nrow(out))
out$npop=rep(npop,nrow(out))

write.table(out,file="SwAsp_94samples.SnIPRE.out.txt",sep="\t",row.name=F,col.name=T,quote=F)







