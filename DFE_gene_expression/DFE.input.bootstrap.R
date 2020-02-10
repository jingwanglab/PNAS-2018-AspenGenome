#! /usr/bin/Rscript --no-save --no-restore

##########################################################################################################
###Main aim: the main aim of this script is to create the input file used in the aspen genome paper
#####first, there should be five gene-character class (5-folders):
		##1. connectivity (high vs. low)
		##2. expression-level (high vs. low)
		##3. expression-variance (high vs. low)
		##4. eGene vs. non-eGene
		##5. core vs. non-core

#####second, within each gene character class, using 4-fold synonymous sites as the neutral sites, I will test the patterns of negative and positive selection in different annotation types:
		##1. 0-fold
		##2. 5-UTR
		##3. 3-UTR
		##4. Intorn
		##5. Upstream
		##6. Downstream

#####third, for annotation type within each gene class, I will build 200 bootstraps to make the plot, and the bootstrap is based on levels, which is different from previous DFE tests when using P. trichocarpa as the reference, previously I used sites as the base to build the bootstrap

library(data.table)
library(dplyr)
library(reshape2)


args=(commandArgs(TRUE))
anno=args[1]   ##annotation categories:0_fold,4_fold,utr3,utr5,intron,upstream,downstram
bootstrap=args[2]    ##bootstrap number 1,2,3,...,200

##########################################################################################################
#####Step1: read in tables published in Plos Genetics, which includes all the information for all expressed genes in the SwAsp samples
eqtls=fread("/proj/b2010014/nobackup/population_genetics/DFE_expression/plos_genetics/S2_file_network_eqtl_gene_statistics.tsv", header=T)  ##the gene-network file
eqtls_all=fread("/proj/b2010014/nobackup/population_genetics/DFE_expression/plos_genetics/S2_file_all_gene_statistics.tsv",header=T)  ##the file is for all gene, with original gene expression and variance without standandization

eqtls_new=eqtls_all[which(eqtls_all$gene %in% eqtls$gene),]
eqtls$orig_expr_mean=eqtls_new$expr_mean
eqtls$orig_expr_var=eqtls_new$expr_var

eqtls_final=eqtls[,c(1,2,9,46,47,15),with=F]  ##this is the final file which I used, which contains the five variables the same as the Figure 6 in the Plos Genetics paper, and will also be used in the aspen genome paper


##########################################################################################################
#####Step2: Divide genes into two classes for different gene classes
###This part is a bit different from the "real" part, since I need to make bootstrap of genes for downstream analysis

#############first define the bootstrap function

f=function(gene){    ###doing bootstrap based on gene levels
        n=c(1:nrow(gene))  ###number of genes in the specific group
        n_boot=sample(n,size=length(n),replace=T)
        boot_gene=as.data.frame(gene[n_boot,])
	colnames(boot_gene)=c("gene")
	return(boot_gene)
}


#############read the dataset

##2.1 connectivity 
##oringal gene files
low_connect_o=filter(eqtls_final,kTotal <= quantile(eqtls_final$kTotal,0.5)) %>% select(gene)
##bootstrap gene file, the same as below
low_connect=f(low_connect_o)

high_connect_o=filter(eqtls_final,kTotal > quantile(eqtls_final$kTotal,0.5)) %>% select(gene)
high_connect=f(high_connect_o)

##2.2 gene expression level
low_exp_o=filter(eqtls_final,orig_expr_mean <= quantile(eqtls_final$orig_expr_mean,0.5)) %>% select(gene)
low_exp=f(low_exp_o)

high_exp_o=filter(eqtls_final,orig_expr_mean > quantile(eqtls_final$orig_expr_mean,0.5)) %>% select(gene)
high_exp=f(high_exp_o)

##2.3 gene expression variance
low_exp_v_o=filter(eqtls_final,orig_expr_var <= quantile(eqtls_final$orig_expr_var,0.5)) %>% select(gene)
low_exp_v=f(low_exp_v_o)

high_exp_v_o=filter(eqtls_final,orig_expr_var > quantile(eqtls_final$orig_expr_var,0.5)) %>% select(gene)
high_exp_v=f(high_exp_v_o)

##2.4 eGene vs. non-eGene
egene_o=filter(eqtls_final,is_egene == "TRUE") %>% select(gene)
egene=f(egene_o)

non_egene_o=filter(eqtls_final,is_egene == "FALSE") %>% select(gene)
non_egene=f(non_egene_o)

##2.5 core vs. non-core
core_o=filter(eqtls_final,is_core_gene == "TRUE") %>% select(gene)
core=f(core_o)

non_core_o=filter(eqtls_final,is_core_gene == "FALSE") %>% select(gene)
non_core=f(non_core_o)

##########################################################################################################
#####Step3: read in the annotation class files, (1) the minor allele count file; (2) the total number of sites file
##3.1 read in the minor allele count file
###always read in the 4-fold file since 4-fold synonymous sites are used as the neutral sites

four_fold_maf=fread("/proj/b2010014/nobackup/population_genetics/DFE_expression/vcf_count_94samples/anno/SwAsp_94samples.polymorphic.maf.4_fold.bed",header=F)
four_fold_maf_unique=unique(four_fold_maf)

##read in the annotation file, which could be 0-fold,utr3,utr5,upstream,downstream
anno_maf=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/vcf_count_94samples/anno/SwAsp_94samples.polymorphic.maf.",anno,".bed",sep="")
anno_maf_file=fread(anno_maf,header=F)
anno_maf_unique=unique(anno_maf_file)


##3.2 read in the total count file
four_fold_total=fread("/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/anno/filter/merge/tremula.4_fold.gene.merge.bed",header=F)
four_fold_n=sum(as.numeric(four_fold_total$V3))-sum(as.numeric(four_fold_total$V2))

####annotated one
anno_total=fread(paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/anno/filter/merge/tremula.",anno,".gene.merge.bed",sep=""),header=F)
anno_n=sum(as.numeric(anno_total$V3))-sum(as.numeric(anno_total$V2))


##########################################################################################################
###Step4 create the folded site freuqency spectrum for different annotation groups
##4.1 connectivity
##4.1.1 4-fold
###count the number of non-polymorphic sites of 4-fold synonymous for low connectivity genes 
##########low connectivity######################
four_fold_total_low_connect=four_fold_total[which(four_fold_total$V4 %in% low_connect$gene),]
four_fold_low_connect_non_poly=sum(as.numeric(four_fold_total_low_connect$V3))-sum(as.numeric(four_fold_total_low_connect$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% low_connect$gene),])
###number of folded minor allele account
four_fold_low_connect_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% low_connect$gene),]$V4)
##build the folded sfs
four_fold_low_connect_sfs=rep(0,189)
four_fold_low_connect_sfs[1]=four_fold_low_connect_non_poly
for (i in 2:95) {
j=i-1
four_fold_low_connect_sfs[i]=four_fold_low_connect_count_table[[j]]}


##########high connectivity######################
four_fold_total_high_connect=four_fold_total[which(four_fold_total$V4 %in% high_connect$gene),]
four_fold_high_connect_non_poly=sum(as.numeric(four_fold_total_high_connect$V3))-sum(as.numeric(four_fold_total_high_connect$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% high_connect$gene),])
###number of folded minor allele account
four_fold_high_connect_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% high_connect$gene),]$V4)
##build the folded sfs
four_fold_high_connect_sfs=rep(0,189)
four_fold_high_connect_sfs[1]=four_fold_high_connect_non_poly
for (i in 2:95) {
j=i-1
four_fold_high_connect_sfs[i]=four_fold_high_connect_count_table[[j]]}


##4.1.2 annotation

##########low connectivity######################
anno_total_low_connect=anno_total[which(anno_total$V4 %in% low_connect$gene),]
anno_low_connect_non_poly=sum(as.numeric(anno_total_low_connect$V3))-sum(as.numeric(anno_total_low_connect$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% low_connect$gene),])
###number of folded minor allele account
anno_low_connect_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% low_connect$gene),]$V4)
##build the folded sfs
anno_low_connect_sfs=rep(0,189)
anno_low_connect_sfs[1]=anno_low_connect_non_poly
for (i in 2:95) {
j=i-1
anno_low_connect_sfs[i]=anno_low_connect_count_table[[j]]}

################################################
##write.table

Outfile_anno_low=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/connectivity/low/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".connectivity.low.real.txt",sep="")

cat("1",file=Outfile_anno_low)
cat("\n",file=Outfile_anno_low,append=T)
cat("188",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(anno_low_connect_sfs,sep=" ",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(four_fold_low_connect_sfs,sep=" ",file=Outfile_anno_low,append=T)


##########high connectivity######################
anno_total_high_connect=anno_total[which(anno_total$V4 %in% high_connect$gene),]
anno_high_connect_non_poly=sum(as.numeric(anno_total_high_connect$V3))-sum(as.numeric(anno_total_high_connect$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% high_connect$gene),])
###number of folded minor allele account
anno_high_connect_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% high_connect$gene),]$V4)
##build the folded sfs
anno_high_connect_sfs=rep(0,189)
anno_high_connect_sfs[1]=anno_high_connect_non_poly
for (i in 2:95) {
j=i-1
anno_high_connect_sfs[i]=anno_high_connect_count_table[[j]]}

###############writing to the table########
################################################
##write.table

Outfile_anno_high=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/connectivity/high/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".connectivity.high.real.txt",sep="")

cat("1",file=Outfile_anno_high)
cat("\n",file=Outfile_anno_high,append=T)
cat("188",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(anno_high_connect_sfs,sep=" ",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(four_fold_high_connect_sfs,sep=" ",file=Outfile_anno_high,append=T)


##########################################################################################################
##4.2 gene expression level

#4.2.1
##########low expression######################
four_fold_total_low_exp=four_fold_total[which(four_fold_total$V4 %in% low_exp$gene),]
four_fold_low_exp_non_poly=sum(as.numeric(four_fold_total_low_exp$V3))-sum(as.numeric(four_fold_total_low_exp$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% low_exp$gene),])
###number of folded minor allele account
four_fold_low_exp_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% low_exp$gene),]$V4)
##build the folded sfs
four_fold_low_exp_sfs=rep(0,189)
four_fold_low_exp_sfs[1]=four_fold_low_exp_non_poly
for (i in 2:95) {
j=i-1
four_fold_low_exp_sfs[i]=four_fold_low_exp_count_table[[j]]}


##########high expression######################
four_fold_total_high_exp=four_fold_total[which(four_fold_total$V4 %in% high_exp$gene),]
four_fold_high_exp_non_poly=sum(as.numeric(four_fold_total_high_exp$V3))-sum(as.numeric(four_fold_total_high_exp$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% high_exp$gene),])
###number of folded minor allele account
four_fold_high_exp_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% high_exp$gene),]$V4)
##build the folded sfs
four_fold_high_exp_sfs=rep(0,189)
four_fold_high_exp_sfs[1]=four_fold_high_exp_non_poly
for (i in 2:95) {
j=i-1
four_fold_high_exp_sfs[i]=four_fold_high_exp_count_table[[j]]}


##4.2.2 annotation

##########low expression######################
anno_total_low_exp=anno_total[which(anno_total$V4 %in% low_exp$gene),]
anno_low_exp_non_poly=sum(as.numeric(anno_total_low_exp$V3))-sum(as.numeric(anno_total_low_exp$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% low_exp$gene),])
###number of folded minor allele account
anno_low_exp_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% low_exp$gene),]$V4)
##build the folded sfs
anno_low_exp_sfs=rep(0,189)
anno_low_exp_sfs[1]=anno_low_exp_non_poly
for (i in 2:95) {
j=i-1
anno_low_exp_sfs[i]=anno_low_exp_count_table[[j]]}

################################################
##write.table

Outfile_anno_low=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_level/low/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".exp_level.low.real.txt",sep="")

cat("1",file=Outfile_anno_low)
cat("\n",file=Outfile_anno_low,append=T)
cat("188",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(anno_low_exp_sfs,sep=" ",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(four_fold_low_exp_sfs,sep=" ",file=Outfile_anno_low,append=T)


##########high expression######################
anno_total_high_exp=anno_total[which(anno_total$V4 %in% high_exp$gene),]
anno_high_exp_non_poly=sum(as.numeric(anno_total_high_exp$V3))-sum(as.numeric(anno_total_high_exp$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% high_exp$gene),])
###number of folded minor allele account
anno_high_exp_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% high_exp$gene),]$V4)
##build the folded sfs
anno_high_exp_sfs=rep(0,189)
anno_high_exp_sfs[1]=anno_high_exp_non_poly
for (i in 2:95) {
j=i-1
anno_high_exp_sfs[i]=anno_high_exp_count_table[[j]]}

###############writing to the table########
################################################
##write.table

Outfile_anno_high=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_level/high/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".exp_level.high.real.txt",sep="")

cat("1",file=Outfile_anno_high)
cat("\n",file=Outfile_anno_high,append=T)
cat("188",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(anno_high_exp_sfs,sep=" ",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(four_fold_high_exp_sfs,sep=" ",file=Outfile_anno_high,append=T)


############################################################################################
##4.3 gene expression variance

#4.3.1
##########low expression######################
four_fold_total_low_exp_v=four_fold_total[which(four_fold_total$V4 %in% low_exp_v$gene),]
four_fold_low_exp_v_non_poly=sum(as.numeric(four_fold_total_low_exp_v$V3))-sum(as.numeric(four_fold_total_low_exp_v$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% low_exp_v$gene),])
###number of folded minor allele account
four_fold_low_exp_v_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% low_exp_v$gene),]$V4)
##build the folded sfs
four_fold_low_exp_v_sfs=rep(0,189)
four_fold_low_exp_v_sfs[1]=four_fold_low_exp_v_non_poly
for (i in 2:95) {
j=i-1
four_fold_low_exp_v_sfs[i]=four_fold_low_exp_v_count_table[[j]]}


##########high expression######################
four_fold_total_high_exp_v=four_fold_total[which(four_fold_total$V4 %in% high_exp_v$gene),]
four_fold_high_exp_v_non_poly=sum(as.numeric(four_fold_total_high_exp_v$V3))-sum(as.numeric(four_fold_total_high_exp_v$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% high_exp_v$gene),])
###number of folded minor allele account
four_fold_high_exp_v_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% high_exp_v$gene),]$V4)
##build the folded sfs
four_fold_high_exp_v_sfs=rep(0,189)
four_fold_high_exp_v_sfs[1]=four_fold_high_exp_v_non_poly
for (i in 2:95) {
j=i-1
four_fold_high_exp_v_sfs[i]=four_fold_high_exp_v_count_table[[j]]}


##4.2.2 annotation

##########low expression######################
anno_total_low_exp_v=anno_total[which(anno_total$V4 %in% low_exp_v$gene),]
anno_low_exp_v_non_poly=sum(as.numeric(anno_total_low_exp_v$V3))-sum(as.numeric(anno_total_low_exp_v$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% low_exp_v$gene),])
###number of folded minor allele account
anno_low_exp_v_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% low_exp_v$gene),]$V4)
##build the folded sfs
anno_low_exp_v_sfs=rep(0,189)
anno_low_exp_v_sfs[1]=anno_low_exp_v_non_poly
for (i in 2:95) {
j=i-1
anno_low_exp_v_sfs[i]=anno_low_exp_v_count_table[[j]]}

################################################
##write.table

Outfile_anno_low=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_variance/low/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".exp_variance.low.real.txt",sep="")

cat("1",file=Outfile_anno_low)
cat("\n",file=Outfile_anno_low,append=T)
cat("188",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(anno_low_exp_v_sfs,sep=" ",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(four_fold_low_exp_v_sfs,sep=" ",file=Outfile_anno_low,append=T)


##########high expression######################
anno_total_high_exp_v=anno_total[which(anno_total$V4 %in% high_exp_v$gene),]
anno_high_exp_v_non_poly=sum(as.numeric(anno_total_high_exp_v$V3))-sum(as.numeric(anno_total_high_exp_v$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% high_exp_v$gene),])
###number of folded minor allele account
anno_high_exp_v_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% high_exp_v$gene),]$V4)
##build the folded sfs
anno_high_exp_v_sfs=rep(0,189)
anno_high_exp_v_sfs[1]=anno_high_exp_v_non_poly
for (i in 2:95) {
j=i-1
anno_high_exp_v_sfs[i]=anno_high_exp_v_count_table[[j]]}

###############writing to the table########
################################################
##write.table

Outfile_anno_high=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_variance/high/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".exp_variance.high.real.txt",sep="")

cat("1",file=Outfile_anno_high)
cat("\n",file=Outfile_anno_high,append=T)
cat("188",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(anno_high_exp_v_sfs,sep=" ",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(four_fold_high_exp_v_sfs,sep=" ",file=Outfile_anno_high,append=T)


#################################################################################################

##4.4 eGene vs. non-eGene

#4.2.1
##########non-egene######################
four_fold_total_non_egene=four_fold_total[which(four_fold_total$V4 %in% non_egene$gene),]
four_fold_non_egene_non_poly=sum(as.numeric(four_fold_total_non_egene$V3))-sum(as.numeric(four_fold_total_non_egene$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% non_egene$gene),])
###number of folded minor allele account
four_fold_non_egene_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% non_egene$gene),]$V4)
##build the folded sfs
four_fold_non_egene_sfs=rep(0,189)
four_fold_non_egene_sfs[1]=four_fold_non_egene_non_poly
for (i in 2:95) {
j=i-1
four_fold_non_egene_sfs[i]=four_fold_non_egene_count_table[[j]]}


##############eGene######################
four_fold_total_egene=four_fold_total[which(four_fold_total$V4 %in% egene$gene),]
four_fold_egene_non_poly=sum(as.numeric(four_fold_total_egene$V3))-sum(as.numeric(four_fold_total_egene$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% egene$gene),])
###number of folded minor allele account
four_fold_egene_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% egene$gene),]$V4)
##build the folded sfs
four_fold_egene_sfs=rep(0,189)
four_fold_egene_sfs[1]=four_fold_egene_non_poly
for (i in 2:95) {
j=i-1
four_fold_egene_sfs[i]=four_fold_egene_count_table[[j]]}

##4.4.2 annotation

##########non_egene##########################
anno_total_non_egene=anno_total[which(anno_total$V4 %in% non_egene$gene),]
anno_non_egene_non_poly=sum(as.numeric(anno_total_non_egene$V3))-sum(as.numeric(anno_total_non_egene$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% non_egene$gene),])
###number of folded minor allele account
anno_non_egene_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% non_egene$gene),]$V4)
##build the folded sfs
anno_non_egene_sfs=rep(0,189)
anno_non_egene_sfs[1]=anno_non_egene_non_poly
for (i in 2:95) {
j=i-1
anno_non_egene_sfs[i]=anno_non_egene_count_table[[j]]}

################################################
##write.table

Outfile_anno_low=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/eGene/low/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".eGene.low.real.txt",sep="")

cat("1",file=Outfile_anno_low)
cat("\n",file=Outfile_anno_low,append=T)
cat("188",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(anno_non_egene_sfs,sep=" ",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(four_fold_non_egene_sfs,sep=" ",file=Outfile_anno_low,append=T)


##########egene##################################
anno_total_egene=anno_total[which(anno_total$V4 %in% egene$gene),]
anno_egene_non_poly=sum(as.numeric(anno_total_egene$V3))-sum(as.numeric(anno_total_egene$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% egene$gene),])
###number of folded minor allele account
anno_egene_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% egene$gene),]$V4)
##build the folded sfs
anno_egene_sfs=rep(0,189)
anno_egene_sfs[1]=anno_egene_non_poly
for (i in 2:95) {
j=i-1
anno_egene_sfs[i]=anno_egene_count_table[[j]]}

###############writing to the table########
################################################
##write.table

Outfile_anno_high=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/eGene/high/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".eGene.high.real.txt",sep="")

cat("1",file=Outfile_anno_high)
cat("\n",file=Outfile_anno_high,append=T)
cat("188",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(anno_egene_sfs,sep=" ",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(four_fold_egene_sfs,sep=" ",file=Outfile_anno_high,append=T)


##4.5 core vs. non-core

#4.5.1
##########non_core######################
four_fold_total_non_core=four_fold_total[which(four_fold_total$V4 %in% non_core$gene),]
four_fold_non_core_non_poly=sum(as.numeric(four_fold_total_non_core$V3))-sum(as.numeric(four_fold_total_non_core$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% non_core$gene),])
###number of folded minor allele account
four_fold_non_core_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% non_core$gene),]$V4)
##build the folded sfs
four_fold_non_core_sfs=rep(0,189)
four_fold_non_core_sfs[1]=four_fold_non_core_non_poly
for (i in 2:95) {
j=i-1
four_fold_non_core_sfs[i]=four_fold_non_core_count_table[[j]]}


##############core######################
four_fold_total_core=four_fold_total[which(four_fold_total$V4 %in% core$gene),]
four_fold_core_non_poly=sum(as.numeric(four_fold_total_core$V3))-sum(as.numeric(four_fold_total_core$V2))-nrow(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% core$gene),])
###number of folded minor allele account
four_fold_core_count_table=table(four_fold_maf_unique[which(four_fold_maf_unique$V5 %in% core$gene),]$V4)
##build the folded sfs
four_fold_core_sfs=rep(0,189)
four_fold_core_sfs[1]=four_fold_core_non_poly
for (i in 2:95) {
j=i-1
four_fold_core_sfs[i]=four_fold_core_count_table[[j]]}

##4.5.2 annotation

##########non_core##########################
anno_total_non_core=anno_total[which(anno_total$V4 %in% non_core$gene),]
anno_non_core_non_poly=sum(as.numeric(anno_total_non_core$V3))-sum(as.numeric(anno_total_non_core$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% non_core$gene),])
###number of folded minor allele account
anno_non_core_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% non_core$gene),]$V4)
##build the folded sfs
anno_non_core_sfs=rep(0,189)
anno_non_core_sfs[1]=anno_non_core_non_poly
for (i in 2:95) {
j=i-1
anno_non_core_sfs[i]=anno_non_core_count_table[[j]]}

################################################
##write.table

Outfile_anno_low=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/core_gene/low/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".core_gene.low.real.txt",sep="")

cat("1",file=Outfile_anno_low)
cat("\n",file=Outfile_anno_low,append=T)
cat("188",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(anno_non_core_sfs,sep=" ",file=Outfile_anno_low,append=T)
cat("\n",file=Outfile_anno_low,append=T)
cat(four_fold_non_core_sfs,sep=" ",file=Outfile_anno_low,append=T)


##########core##################################
anno_total_core=anno_total[which(anno_total$V4 %in% core$gene),]
anno_core_non_poly=sum(as.numeric(anno_total_core$V3))-sum(as.numeric(anno_total_core$V2))-nrow(anno_maf_unique[which(anno_maf_unique$V5 %in% core$gene),])
###number of folded minor allele account
anno_core_count_table=table(anno_maf_unique[which(anno_maf_unique$V5 %in% core$gene),]$V4)
##build the folded sfs
anno_core_sfs=rep(0,189)
anno_core_sfs[1]=anno_core_non_poly
for (i in 2:95) {
j=i-1
anno_core_sfs[i]=anno_core_count_table[[j]]}

###############writing to the table########
################################################
##write.table

Outfile_anno_high=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/core_gene/high/",anno,"/bootstrap/bootstrap",bootstrap,"/",anno,".core_gene.high.real.txt",sep="")

cat("1",file=Outfile_anno_high)
cat("\n",file=Outfile_anno_high,append=T)
cat("188",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(anno_core_sfs,sep=" ",file=Outfile_anno_high,append=T)
cat("\n",file=Outfile_anno_high,append=T)
cat(four_fold_core_sfs,sep=" ",file=Outfile_anno_high,append=T)



