#! /usr/bin/Rscript --no-save --no-restore

##########################################################################################################
###Main aim: the main aim of this script is to create the input file used in the aspen genome paper
##This script is to create the input file for omega from DFE to run

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

#####third, for annotation type within each gene class, I will build 200 bootstraps to make the plot (this part maybe included in other scripts...)

library(data.table)
library(dplyr)


args=(commandArgs(TRUE))
anno=args[1]

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

##2.1 connectivity 
low_connectivity=filter(eqtls_final,kTotal <= quantile(eqtls_final$kTotal,0.5)) %>% select(gene)
high_connectivity=filter(eqtls_final,kTotal > quantile(eqtls_final$kTotal,0.5)) %>% select(gene)

##2.2 gene expression level
low_exp_level=filter(eqtls_final,orig_expr_mean <= quantile(eqtls_final$orig_expr_mean,0.5)) %>% select(gene)
high_exp_level=filter(eqtls_final,orig_expr_mean > quantile(eqtls_final$orig_expr_mean,0.5)) %>% select(gene)

##2.3 gene expression variance
low_exp_variance=filter(eqtls_final,orig_expr_var <= quantile(eqtls_final$orig_expr_var,0.5)) %>% select(gene)
high_exp_variance=filter(eqtls_final,orig_expr_var > quantile(eqtls_final$orig_expr_var,0.5)) %>% select(gene)

##2.4 eGene vs. non-eGene
high_eGene=filter(eqtls_final,is_egene == "TRUE") %>% select(gene)
low_eGene=filter(eqtls_final,is_egene == "FALSE") %>% select(gene)

##2.5 core vs. non-core
high_core_gene=filter(eqtls_final,is_core_gene == "TRUE") %>% select(gene)
low_core_gene=filter(eqtls_final,is_core_gene == "FALSE") %>% select(gene)


##########################################################################################################
#####Step3: read in the annotation class files, (1) the trichocarpa genotype file
##3.1 read in the trichocarpa file, which contains the sites where the sample of trichocarpa is the same as tremula, the one which is heterzygotes, and the one that is homozygotes different from tremula
###always read in the 4-fold file since 4-fold synonymous sites are used as the neutral sites

four_fold_tr=fread("/proj/b2010014/nobackup/population_genetics/DFE_expression/outgroup_vcf/trichocarpa/bed_intersect/trichocarpa.frq.filter.4_fold.bed",header=F)

##read in the annotation file, which could be 0-fold,utr3,utr5,upstream,downstream
anno_tr_file=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/outgroup_vcf/trichocarpa/bed_intersect/trichocarpa.frq.filter.",anno,".bed",sep="")
anno_tr=fread(anno_tr_file,header=F)

##########################################################################################################
###Step4 create the divergence file for different annotation groups
##4.1 connectivity

##4.1.1 4-fold
###count the number of all mapped sites of 4-fold synonymou and the number of fixed sites between the individual of trichocarpar and the reference genome of tremulas for low connectivity genes
##########low connectivity######################
##four fold sites used as the neutral sites
four_fold_tr_low_connectivity=four_fold_tr[which(four_fold_tr$V6 %in% low_connectivity$gene),]
four_fold_tr_low_connectivity_t=nrow(four_fold_tr_low_connectivity)   ##total number of mapped sites for four fold in low connectivity genes
four_fold_tr_low_connectivity_d=length(which(four_fold_tr_low_connectivity$V4==2 & four_fold_tr_low_connectivity$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

#4.1.2 anno
##all other functional sites 
anno_tr_low_connectivity=anno_tr[which(anno_tr$V6 %in% low_connectivity$gene),]
anno_tr_low_connectivity_t=nrow(anno_tr_low_connectivity)   ##total number of mapped sites for four fold in low connectivity genes
anno_tr_low_connectivity_d=length(which(anno_tr_low_connectivity$V4==2 & anno_tr_low_connectivity$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_low_connectivity=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/connectivity/low/",anno,"/real/",anno,".connectivity.low.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_low_connectivity_t,anno_tr_low_connectivity_d,file=Outfile_low_connectivity,sep=" ")
cat("\n",file=Outfile_low_connectivity,append=T)
cat("0",four_fold_tr_low_connectivity_t,four_fold_tr_low_connectivity_d,file=Outfile_low_connectivity,sep=" ",append=T)

##########################################################################################
##########high connectivity######################
##four fold sites used as the neutral sites
four_fold_tr_high_connectivity=four_fold_tr[which(four_fold_tr$V6 %in% high_connectivity$gene),]
four_fold_tr_high_connectivity_t=nrow(four_fold_tr_high_connectivity)   ##total number of mapped sites for four fold in high connectivity genes
four_fold_tr_high_connectivity_d=length(which(four_fold_tr_high_connectivity$V4==2 & four_fold_tr_high_connectivity$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

##all other functional sites
anno_tr_high_connectivity=anno_tr[which(anno_tr$V6 %in% high_connectivity$gene),]
anno_tr_high_connectivity_t=nrow(anno_tr_high_connectivity)   ##total number of mapped sites for four fold in high connectivity genes
anno_tr_high_connectivity_d=length(which(anno_tr_high_connectivity$V4==2 & anno_tr_high_connectivity$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_high_connectivity=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/connectivity/high/",anno,"/real/",anno,".connectivity.high.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_high_connectivity_t,anno_tr_high_connectivity_d,file=Outfile_high_connectivity,sep=" ")
cat("\n",file=Outfile_high_connectivity,append=T)
cat("0",four_fold_tr_high_connectivity_t,four_fold_tr_high_connectivity_d,file=Outfile_high_connectivity,sep=" ",append=T)


##########################################################################################
##4.2 exp_level

##4.1.1 4-fold
###count the number of all mapped sites of 4-fold synonymou and the number of fixed sites between the individual of trichocarpar and the reference genome of tremulas for low exp_level genes
##########low exp_level######################
##four fold sites used as the neutral sites
four_fold_tr_low_exp_level=four_fold_tr[which(four_fold_tr$V6 %in% low_exp_level$gene),]
four_fold_tr_low_exp_level_t=nrow(four_fold_tr_low_exp_level)   ##total number of mapped sites for four fold in low exp_level genes
four_fold_tr_low_exp_level_d=length(which(four_fold_tr_low_exp_level$V4==2 & four_fold_tr_low_exp_level$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

#4.1.2 anno
##all other functional sites
anno_tr_low_exp_level=anno_tr[which(anno_tr$V6 %in% low_exp_level$gene),]
anno_tr_low_exp_level_t=nrow(anno_tr_low_exp_level)   ##total number of mapped sites for four fold in low exp_level genes
anno_tr_low_exp_level_d=length(which(anno_tr_low_exp_level$V4==2 & anno_tr_low_exp_level$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_low_exp_level=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_level/low/",anno,"/real/",anno,".exp_level.low.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_low_exp_level_t,anno_tr_low_exp_level_d,file=Outfile_low_exp_level,sep=" ")
cat("\n",file=Outfile_low_exp_level,append=T)
cat("0",four_fold_tr_low_exp_level_t,four_fold_tr_low_exp_level_d,file=Outfile_low_exp_level,sep=" ",append=T)

##########################################################################################
##########high exp_level######################
##four fold sites used as the neutral sites
four_fold_tr_high_exp_level=four_fold_tr[which(four_fold_tr$V6 %in% high_exp_level$gene),]
four_fold_tr_high_exp_level_t=nrow(four_fold_tr_high_exp_level)   ##total number of mapped sites for four fold in high exp_level genes
four_fold_tr_high_exp_level_d=length(which(four_fold_tr_high_exp_level$V4==2 & four_fold_tr_high_exp_level$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

##all other functional sites
anno_tr_high_exp_level=anno_tr[which(anno_tr$V6 %in% high_exp_level$gene),]
anno_tr_high_exp_level_t=nrow(anno_tr_high_exp_level)   ##total number of mapped sites for four fold in high exp_level genes
anno_tr_high_exp_level_d=length(which(anno_tr_high_exp_level$V4==2 & anno_tr_high_exp_level$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_high_exp_level=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_level/high/",anno,"/real/",anno,".exp_level.high.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_high_exp_level_t,anno_tr_high_exp_level_d,file=Outfile_high_exp_level,sep=" ")
cat("\n",file=Outfile_high_exp_level,append=T)
cat("0",four_fold_tr_high_exp_level_t,four_fold_tr_high_exp_level_d,file=Outfile_high_exp_level,sep=" ",append=T)


#####################################################################################################################
##4.1 exp_variance

##4.1.1 4-fold
###count the number of all mapped sites of 4-fold synonymou and the number of fixed sites between the individual of trichocarpar and the reference genome of tremulas for low exp_variance genes
##########low exp_variance######################
##four fold sites used as the neutral sites
four_fold_tr_low_exp_variance=four_fold_tr[which(four_fold_tr$V6 %in% low_exp_variance$gene),]
four_fold_tr_low_exp_variance_t=nrow(four_fold_tr_low_exp_variance)   ##total number of mapped sites for four fold in low exp_variance genes
four_fold_tr_low_exp_variance_d=length(which(four_fold_tr_low_exp_variance$V4==2 & four_fold_tr_low_exp_variance$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

#4.1.2 anno
##all other functional sites
anno_tr_low_exp_variance=anno_tr[which(anno_tr$V6 %in% low_exp_variance$gene),]
anno_tr_low_exp_variance_t=nrow(anno_tr_low_exp_variance)   ##total number of mapped sites for four fold in low exp_variance genes
anno_tr_low_exp_variance_d=length(which(anno_tr_low_exp_variance$V4==2 & anno_tr_low_exp_variance$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_low_exp_variance=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_variance/low/",anno,"/real/",anno,".exp_variance.low.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_low_exp_variance_t,anno_tr_low_exp_variance_d,file=Outfile_low_exp_variance,sep=" ")
cat("\n",file=Outfile_low_exp_variance,append=T)
cat("0",four_fold_tr_low_exp_variance_t,four_fold_tr_low_exp_variance_d,file=Outfile_low_exp_variance,sep=" ",append=T)


##########################################################################################
##########high exp_variance######################
##four fold sites used as the neutral sites
four_fold_tr_high_exp_variance=four_fold_tr[which(four_fold_tr$V6 %in% high_exp_variance$gene),]
four_fold_tr_high_exp_variance_t=nrow(four_fold_tr_high_exp_variance)   ##total number of mapped sites for four fold in high exp_variance genes
four_fold_tr_high_exp_variance_d=length(which(four_fold_tr_high_exp_variance$V4==2 & four_fold_tr_high_exp_variance$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

##all other functional sites
anno_tr_high_exp_variance=anno_tr[which(anno_tr$V6 %in% high_exp_variance$gene),]
anno_tr_high_exp_variance_t=nrow(anno_tr_high_exp_variance)   ##total number of mapped sites for four fold in high exp_variance genes
anno_tr_high_exp_variance_d=length(which(anno_tr_high_exp_variance$V4==2 & anno_tr_high_exp_variance$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_high_exp_variance=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/exp_variance/high/",anno,"/real/",anno,".exp_variance.high.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_high_exp_variance_t,anno_tr_high_exp_variance_d,file=Outfile_high_exp_variance,sep=" ")
cat("\n",file=Outfile_high_exp_variance,append=T)
cat("0",four_fold_tr_high_exp_variance_t,four_fold_tr_high_exp_variance_d,file=Outfile_high_exp_variance,sep=" ",append=T)


###############################################################################################################
##4.4 eGene

##4.1.1 4-fold
###count the number of all mapped sites of 4-fold synonymou and the number of fixed sites between the individual of trichocarpar and the reference genome of tremulas for low eGene genes
##########low eGene######################
##four fold sites used as the neutral sites
four_fold_tr_low_eGene=four_fold_tr[which(four_fold_tr$V6 %in% low_eGene$gene),]
four_fold_tr_low_eGene_t=nrow(four_fold_tr_low_eGene)   ##total number of mapped sites for four fold in low eGene genes
four_fold_tr_low_eGene_d=length(which(four_fold_tr_low_eGene$V4==2 & four_fold_tr_low_eGene$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

#4.1.2 anno
##all other functional sites
anno_tr_low_eGene=anno_tr[which(anno_tr$V6 %in% low_eGene$gene),]
anno_tr_low_eGene_t=nrow(anno_tr_low_eGene)   ##total number of mapped sites for four fold in low eGene genes
anno_tr_low_eGene_d=length(which(anno_tr_low_eGene$V4==2 & anno_tr_low_eGene$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_low_eGene=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/eGene/low/",anno,"/real/",anno,".eGene.low.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_low_eGene_t,anno_tr_low_eGene_d,file=Outfile_low_eGene,sep=" ")
cat("\n",file=Outfile_low_eGene,append=T)
cat("0",four_fold_tr_low_eGene_t,four_fold_tr_low_eGene_d,file=Outfile_low_eGene,sep=" ",append=T)


##########################################################################################
##########high eGene######################
##four fold sites used as the neutral sites
four_fold_tr_high_eGene=four_fold_tr[which(four_fold_tr$V6 %in% high_eGene$gene),]
four_fold_tr_high_eGene_t=nrow(four_fold_tr_high_eGene)   ##total number of mapped sites for four fold in high eGene genes
four_fold_tr_high_eGene_d=length(which(four_fold_tr_high_eGene$V4==2 & four_fold_tr_high_eGene$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

##all other functional sites
anno_tr_high_eGene=anno_tr[which(anno_tr$V6 %in% high_eGene$gene),]
anno_tr_high_eGene_t=nrow(anno_tr_high_eGene)   ##total number of mapped sites for four fold in high eGene genes
anno_tr_high_eGene_d=length(which(anno_tr_high_eGene$V4==2 & anno_tr_high_eGene$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_high_eGene=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/eGene/high/",anno,"/real/",anno,".eGene.high.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_high_eGene_t,anno_tr_high_eGene_d,file=Outfile_high_eGene,sep=" ")
cat("\n",file=Outfile_high_eGene,append=T)
cat("0",four_fold_tr_high_eGene_t,four_fold_tr_high_eGene_d,file=Outfile_high_eGene,sep=" ",append=T)


#############################################################################################################
##4.1 core_gene

##4.1.1 4-fold
###count the number of all mapped sites of 4-fold synonymou and the number of fixed sites between the individual of trichocarpar and the reference genome of tremulas for low core_gene genes
##########low core_gene######################
##four fold sites used as the neutral sites
four_fold_tr_low_core_gene=four_fold_tr[which(four_fold_tr$V6 %in% low_core_gene$gene),]
four_fold_tr_low_core_gene_t=nrow(four_fold_tr_low_core_gene)   ##total number of mapped sites for four fold in low core_gene genes
four_fold_tr_low_core_gene_d=length(which(four_fold_tr_low_core_gene$V4==2 & four_fold_tr_low_core_gene$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

#4.1.2 anno
##all other functional sites
anno_tr_low_core_gene=anno_tr[which(anno_tr$V6 %in% low_core_gene$gene),]
anno_tr_low_core_gene_t=nrow(anno_tr_low_core_gene)   ##total number of mapped sites for four fold in low core_gene genes
anno_tr_low_core_gene_d=length(which(anno_tr_low_core_gene$V4==2 & anno_tr_low_core_gene$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_low_core_gene=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/core_gene/low/",anno,"/real/",anno,".core_gene.low.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_low_core_gene_t,anno_tr_low_core_gene_d,file=Outfile_low_core_gene,sep=" ")
cat("\n",file=Outfile_low_core_gene,append=T)
cat("0",four_fold_tr_low_core_gene_t,four_fold_tr_low_core_gene_d,file=Outfile_low_core_gene,sep=" ",append=T)


##########high core_gene######################
##four fold sites used as the neutral sites
four_fold_tr_high_core_gene=four_fold_tr[which(four_fold_tr$V6 %in% high_core_gene$gene),]
four_fold_tr_high_core_gene_t=nrow(four_fold_tr_high_core_gene)   ##total number of mapped sites for four fold in high core_gene genes
four_fold_tr_high_core_gene_d=length(which(four_fold_tr_high_core_gene$V4==2 & four_fold_tr_high_core_gene$V5==0))  ##the total number fixed differences between one sample of trichocarpa and the reference genome of tremula, here V3==2 indicates there are two alleles, and V5==0 indicates there is no allele that were same as the reference genome

##all other functional sites
anno_tr_high_core_gene=anno_tr[which(anno_tr$V6 %in% high_core_gene$gene),]
anno_tr_high_core_gene_t=nrow(anno_tr_high_core_gene)   ##total number of mapped sites for four fold in high core_gene genes
anno_tr_high_core_gene_d=length(which(anno_tr_high_core_gene$V4==2 & anno_tr_high_core_gene$V5==0))  ##the total number fixed differences between one sample of


##########write.table
Outfile_high_core_gene=paste("/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/core_gene/high/",anno,"/real/",anno,".core_gene.high.omega.txt",sep="")

###cat the number into the table
cat("1",anno_tr_high_core_gene_t,anno_tr_high_core_gene_d,file=Outfile_high_core_gene,sep=" ")
cat("\n",file=Outfile_high_core_gene,append=T)
cat("0",four_fold_tr_high_core_gene_t,four_fold_tr_high_core_gene_d,file=Outfile_high_core_gene,sep=" ",append=T)




