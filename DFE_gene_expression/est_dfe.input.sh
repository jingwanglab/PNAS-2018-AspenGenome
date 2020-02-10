#! /bin/bash -l

#SBATCH -A b2011141
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o est_dfe.input.out
#SBATCH -e est_dfe.input.err
#SBATCH -J est_dfe.input.job
#SBATCH -t 24:00:00

######Main aim: the main aim of this script is used to create the input file for est_dfe related softwares to run, the output file is as following for different site categories: 0-fold, 4-fold, 3UTR, 5UTR, Intronic, Upstream, Downstream:
	#Scaffold	Start	End	Minor_allele_count(all counts for 188 haplotypes)	Gene
	#SNPs...
	#Then, the non-polymorphic sites are calculated by substracting the total number of sites from .bed file by the polymorphic sites from vcf file in each annotation group

module load bioinfo-tools
module load BEDTools

step=$1  ##indicates which step to perfrom

###############################################################################################################################################
#####Step0: I just note that some of the annotation category file I used are from Potra_v1.0, so before all downstream steps, I need to first extract all annotation categories from the new annotation Potra v1.1
#but, since the annotation categories (4_fold and 0_fold) files have been created when I work on the co-expression paper that already published in PlosGenetics. So I can just copy those related annotation files to the corresponding folders


if [ "$step" == "0" ]; then
potra_gene_gff="/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/anno/Potra01-gene-representative.gff3"
potra_regulatory_gff="/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/anno/Potra01-1000-regulatory-regions.gff3"

gff_dir=`dirname $potra_gene_gff`

###Step 0.1 extract the annotation categories 

grep "upstream"  $potra_regulatory_gff > $gff_dir/upstream.gff3
grep "downstream"  $potra_regulatory_gff > $gff_dir/downstream.gff3
grep "intron" $potra_gene_gff > $gff_dir/intron.gff3
grep "three_prime_UTR" $potra_gene_gff > $gff_dir/utr3.gff3
grep "five_prime_UTR" $potra_gene_gff > $gff_dir/utr5.gff3
grep "gene" $potra_gene_gff > $gff_dir/gene.gff3
bedtools intersect -a $gff_dir/0_fold.bed -b $gff_dir/gene.gff3 -wb > $gff_dir/0_fold.gene.bed
bedtools intersect -a $gff_dir/4_fold.bed -b $gff_dir/gene.gff3 -wb > $gff_dir/4_fold.gene.bed


###############################################################################################################################################
#####Step1: Using BEDTools to intersect the .bed file after several steps of filtering and the .bed file for each annotation category

elif [ "$step" == "1" ]; then

#anno=$2  ##this is the name of annotation type, e.g. 0_fold, 3UTR, 5UTR, 4_fold, intron, regulatory.downstream, regulatory.upstream

bed_filter="/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/Potra.SwAsp.filter.bed"
bed_dir=`dirname $bed_filter`

#anno_bed="$bed_dir/anno/$anno.bed"

anno_filter_bed="$bed_dir/anno/filter"

if [ ! -d "$anno_filter_bed" ]; then
mkdir -p $anno_filter_bed
fi

#bedtools intersect -a $anno_bed -b $bed_filter > $anno_filter_bed/tremula.$anno.filter.bed
#bedtools intersect -a $bed_filter -b $bed_dir/anno/0_fold.gene.bed -wb > $anno_filter_bed/tremula.0_fold.filter.bed
#bedtools intersect -a $bed_filter -b $bed_dir/anno/4_fold.gene.bed -wb > $anno_filter_bed/tremula.4_fold.filter.bed
#bedtools intersect -a $bed_filter -b $bed_dir/anno/intron.gff3 -wb > $anno_filter_bed/tremula.intron.filter.bed
#bedtools intersect -a $bed_filter -b $bed_dir/anno/utr3.gff3 -wb > $anno_filter_bed/tremula.utr3.filter.bed
#bedtools intersect -a $bed_filter -b $bed_dir/anno/utr5.gff3 -wb > $anno_filter_bed/tremula.utr5.filter.bed
bedtools intersect -a $bed_filter -b $bed_dir/anno/upstream.gff3 -wb > $anno_filter_bed/tremula.upstream.filter.bed
bedtools intersect -a $bed_filter -b $bed_dir/anno/downstream.gff3 -wb > $anno_filter_bed/tremula.downstream.filter.bed


elif [ "$step" == "2" ]; then

##########################################################
###Step2: extract all gene names corresponding to the bed file
###This script is used to only extract the gene name of the bed files
module load BEDTools

bed_filter="/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/Potra.SwAsp.filter.bed"
bed_dir=`dirname $bed_filter`
anno_filter_bed="$bed_dir/anno/filter"

#cut -f 1,2,3,15 $anno_filter_bed/tremula.4_fold.filter.bed |cut -f 1 -d ";" |sed 's/ID=//g' > $anno_filter_bed/tremula.4_fold.gene.bed && rm $anno_filter_bed/tremula.4_fold.filter.bed
#cut -f 1,2,3,15 $anno_filter_bed/tremula.0_fold.filter.bed |cut -f 1 -d ";" |sed 's/ID=//g' > $anno_filter_bed/tremula.0_fold.gene.bed && rm $anno_filter_bed/tremula.0_fold.filter.bed
#cut -f 1,2,3,12 $anno_filter_bed/tremula.downstream.filter.bed |cut -f 1 -d ";" |sed 's/ID=//g' > $anno_filter_bed/tremula.donwstream.gene.bed && rm $anno_filter_bed/tremula.downstream.filter.bed
#cut -f 1,2,3,12 $anno_filter_bed/tremula.upstream.filter.bed |cut -f 1 -d ";" |sed 's/ID=//g' > $anno_filter_bed/tremula.upstream.gene.bed && rm $anno_filter_bed/tremula.upstream.filter.bed
#cut -f 1,2,3,12 $anno_filter_bed/tremula.utr3.filter.bed |cut -f 1 -d ";" |cut -f 1 -d "." |sed 's/ID=//g' > $anno_filter_bed/tremula.utr3.gene.bed && rm $anno_filter_bed/tremula.utr3.filter.bed
#cut -f 1,2,3,12 $anno_filter_bed/tremula.utr5.filter.bed |cut -f 1 -d ";" |cut -f 1 -d "." |sed 's/ID=//g' > $anno_filter_bed/tremula.utr5.gene.bed && rm $anno_filter_bed/tremula.utr5.filter.bed
#cut -f 1,2,3,12 $anno_filter_bed/tremula.intron.filter.bed |cut -f 1 -d "." |sed 's/Parent=//g' > $anno_filter_bed/tremula.intron.gene.bed && rm $anno_filter_bed/tremula.intron.filter.bed

###############################################
anno=$2
anno_merge_bed="$bed_dir/anno/filter/merge"
sort -k1,1 -k2,2n $anno_filter_bed/tremula.$anno.gene.bed > $anno_merge_bed/tremula.$anno.gene.bed 
##if there is overlapped genes, then only keep the first one
bedtools merge -i $anno_merge_bed/tremula.$anno.gene.bed -c 4 -o collapse |cut -f 1 -d "," > $anno_merge_bed/tremula.$anno.gene.merge.bed && rm $anno_merge_bed/tremula.$anno.gene.bed


elif [ "$step" == "3" ]; then
###############################################################################################################################################
#####Step3: transfer the count file to .bed file and then use bedtools intersect to extract the SNP information for each annotation categories

##############################################
#step3.1 first use vcftools --count2 to output the reference and non-reference allele counts across samples, then use perl scripts to transfer the count fileto bed file for downstream intersecting 


count="/proj/b2010014/nobackup/population_genetics/DFE_expression/vcf_count_94samples/SwAsp_94samples.filter.frq.count"  ##count file directly created by vcftools
count_dir=`dirname $count`

count_to_bed="/proj/b2011141/pipeline/perl/genomic_adaptation_paper/aspen_genome_DFE/vcf_count_bed.pl"  ##perl script to transfer count file to .bed file
##use perl to transfer the count to bed file
#perl $count_to_bed $count

##############################################
#step 2.2 use bedtools intersect to extract the SNP Minor allele count information for each annotation categories

anno=$2  ##this is the name of annotation type, e.g. 0_fold, 3UTR, 5UTR, 4_fold, intron, regulatory.downstream, regulatory.upstream

bed_filter="/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/Potra.SwAsp.filter.bed"
bed_dir=`dirname $bed_filter`
anno_filter_bed="$bed_dir/anno/filter/tremula.$anno.gene.bed"

count_anno="$count_dir/anno"

if [ ! -d "$count_anno" ]; then
mkdir -p $count_anno
fi

#bedtools intersect -a $count_dir/SwAsp_94samples.filter.frq.bed -b $anno_filter_bed -wb | cut -f 1,2,3,4,8 > $count_anno/SwAsp_94samples.polymorphic.maf.$anno.bed
#sort -k1,1 -k2,2n $count_anno/SwAsp_94samples.polymorphic.maf.$anno.bed > $count_anno/temp.$anno && mv $count_anno/temp.$anno $count_anno/SwAsp_94samples.polymorphic.maf.$anno.bed


elif [ "$step" == "4" ]; then
###############################################################################################################################################
###Step 4: create the output folders for each gene class and annotation categories
OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out"

anno=$2  ##annotation: 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene 
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for category in {low,high}
do

Out_anno="$OutDir/$gene_class/$category/$anno/real"

if [ ! -d "$Out_anno" ]; then
mkdir -p $Out_anno
fi

done
done


dfe_input="/proj/b2011141/pipeline/R/DFE/aspen_genome_paper/DFE.input.R"

anno=$2  ##annotation: 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

Rscript $dfe_input $anno


###############################################################################################################################################
elif [ "$step" == "5" ]; then

###Step 5: same as step4, just do 200 times of bootstrap here
OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out"

anno=$2  ##annotation: 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for category in {low,high}
do
for bootstrap_n in {1..200}
do

Out_anno="$OutDir/$gene_class/$category/$anno/bootstrap/bootstrap$bootstrap_n"

if [ ! -d "$Out_anno" ]; then
mkdir -p $Out_anno
fi

done
done
done


dfe_bootstrap_input="/proj/b2011141/pipeline/R/DFE/aspen_genome_paper/DFE.input.bootstrap.R"

anno=$2  ##annotation: 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

for boot in {1..200}
do
Rscript $dfe_bootstrap_input $anno $boot
done


###############################################################################################################################################
elif [ "$step" == "6" ]; then
################becauase we have the vcf file that use one individual of P.trichcoarpa to map to the P. tremula reference genome for all sites, I used this vcf file to calculate how many sites that are fixed between P. tremula and P. trichocarpa

#####5.1 first, we use vcf file to output the sites that have relatively good quality,  and have reads mapped to the reference genome, with read depth between 5 and 70 in order to remove those "problem" sites

module load vcftools

trichocarpa_ug_vcf="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/ancestral_alleles/UG_all_SNPS/trichocarpa/SRR1571343.gatk.ug.vcf.gz"

out="/proj/b2010014/nobackup/population_genetics/DFE_expression/outgroup_vcf/trichocarpa"
vcftools --gzvcf $trichocarpa_ug_vcf --remove-indels --max-missing 1 --minQ 30 --minDP 5 --maxDP 70 --counts2 --out $out/trichocarpa

######5.2 use perl script to transfer the output .count file to .bed file for downstream analysis, that use BEDTools to extract the filtered regions for different annotation categories, e.g. 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

count_to_bed="/proj/b2011141/pipeline/DFE/aspen_genome_paper/DFE.count_to_bed.pl"
trichocarpa_count="/proj/b2010014/nobackup/population_genetics/DFE_expression/outgroup_vcf/trichocarpa/trichocarpa.frq.count"

perl $count_to_bed $trichocarpa_count


###############################################################################################################################################
elif [ "$step" == "7" ]; then

###In the step6, use BEDTools to intersect the output .bed file from counts file which includes the information of mapped sites of one sample of P. trichocarpa and also the genotype information of the P. trichocarpa samples
anno=$2  ###the type of annotation categories...

trichocarpa_bed="/proj/b2010014/nobackup/population_genetics/DFE_expression/outgroup_vcf/trichocarpa/trichocarpa.frq.bed"
bed_dir=`dirname $trichocarpa_bed`
Input_bed="/proj/b2010014/nobackup/population_genetics/DFE_expression/bed_94samples/anno/filter/merge"
anno_bed="$Input_bed/tremula.$anno.gene.merge.bed"

OutDir="$bed_dir/bed_intersect"

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

bedtools intersect -a $trichocarpa_bed -b $anno_bed -wa -wb |cut -f 1-5,9 > $OutDir/trichocarpa.frq.filter.$anno.bed


###############################################################################################################################################
elif [ "$step" == "8" ]; then
###Step 8: create the output folders for each gene class and annotation categories for the divergence file
OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out"

anno=$2  ##annotation: 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for category in {low,high}
do

Out_anno="$OutDir/$gene_class/$category/$anno/real"

if [ ! -d "$Out_anno" ]; then
mkdir -p $Out_anno
fi

done
done


dfe_omega_input="/proj/b2011141/pipeline/R/DFE/aspen_genome_paper/DFE.input.omega.R"

Rscript $dfe_omega_input $anno


###############################################################################################################################################
elif [ "$step" == "9" ]; then

###Step 9: same as step8, just do 200 times of bootstrap here
OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out"

anno=$2  ##annotation: 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for category in {low,high}
do
for bootstrap_n in {1..200}
do

Out_anno="$OutDir/$gene_class/$category/$anno/bootstrap/bootstrap$bootstrap_n"

if [ ! -d "$Out_anno" ]; then
mkdir -p $Out_anno
fi

done
done
done


dfe_bootstrap_omega="/proj/b2011141/pipeline/R/DFE/aspen_genome_paper/DFE.input.omega.bootstrap.R"

anno=$2  ##annotation: 0_fold,4_fold,utr3,utr5,intron,upstream,downstream

for boot in {1..200}
do
Rscript $dfe_bootstrap_omega $anno $boot
done


##############################################################################
###Step10: because I found the non-coding regions (especially upstream and downstream) regions of different classes fo genes often show the opposite DFE paper compared to 0-fold nonsynonymous sites,so, I will compare the DFE between two pairs of (high vs. low) for 4_fold synonymous sites of different gene classes to check whether this could just be caused by different selection pattern in synonymous sites (which may not really reflect the true "neutral" sites)

elif [ "$step" == "10" ]; then

OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out"

for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do

######first creating the relevant folders
Out_4_fold_real="$OutDir/$gene_class/4_fold_compare/4_fold/real"

if [ ! -d "$Out_4_fold_real" ]; then
mkdir -p $Out_4_fold_real
fi

for bootstrap_n in {1..200}
do

Out_4_fold_boot="$OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$bootstrap_n"

if [ ! -d "$Out_4_fold_boot" ]; then
mkdir -p $Out_4_fold_boot
fi

done
done

#################copy the relevant file to the corresponding folders 
###create the .real.txt file
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
##4_fold.connectivity.4_fold_compare.real.txt
sed '3d' $OutDir/$gene_class/high/0_fold/real/0_fold.$gene_class.high.real.txt > $OutDir/$gene_class/4_fold_compare/4_fold/real/4_fold.$gene_class.4_fold_compare.real.txt
echo -e "\n" >> $OutDir/$gene_class/4_fold_compare/4_fold/real/4_fold.$gene_class.4_fold_compare.real.txt
tail -n 1 $OutDir/$gene_class/low/0_fold/real/0_fold.$gene_class.low.real.txt >> $OutDir/$gene_class/4_fold_compare/4_fold/real/4_fold.$gene_class.4_fold_compare.real.txt
sed '4d' $OutDir/$gene_class/4_fold_compare/4_fold/real/4_fold.$gene_class.4_fold_compare.real.txt > temp && mv temp $OutDir/$gene_class/4_fold_compare/4_fold/real/4_fold.$gene_class.4_fold_compare.real.txt
done
####################################################################
###create the .omega.txt file
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
##4_fold.connectivity.4_fold_compare.real.txt
sed '1d' $OutDir/$gene_class/high/0_fold/real/0_fold.$gene_class.high.omega.txt | awk '$1=1' > $OutDir/$gene_class/4_fold_compare/4_fold/real/4_fold.$gene_class.4_fold_compare.omega.txt
tail -n 1 $OutDir/$gene_class/low/0_fold/real/0_fold.$gene_class.low.omega.txt  >> $OutDir/$gene_class/4_fold_compare/4_fold/real/4_fold.$gene_class.4_fold_compare.omega.txt

done

####################################################################
###extract the bootstrap results

###create the bootstrap .real.txt file
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
##4_fold.connectivity.4_fold_compare.real.txt
for boot in {1..200}
do
sed '3d' $OutDir/$gene_class/high/0_fold/bootstrap/bootstrap$boot/0_fold.$gene_class.high.real.txt > $OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$boot/4_fold.$gene_class.4_fold_compare.real.txt
echo -e "\n" >> $OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$boot/4_fold.$gene_class.4_fold_compare.real.txt
tail -n 1 $OutDir/$gene_class/low/0_fold/bootstrap/bootstrap$boot/0_fold.$gene_class.low.real.txt >> $OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$boot/4_fold.$gene_class.4_fold_compare.real.txt
sed '4d' $OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$boot/4_fold.$gene_class.4_fold_compare.real.txt > temp && mv temp $OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$boot/4_fold.$gene_class.4_fold_compare.real.txt
done
done
####################################################################
###create the .omega.txt file
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
##4_fold.connectivity.4_fold_compare.real.txt
sed '1d' $OutDir/$gene_class/high/0_fold/bootstrap/bootstrap$boot/0_fold.$gene_class.high.omega.txt | awk '$1=1' > $OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$boot/4_fold.$gene_class.4_fold_compare.omega.txt
tail -n 1 $OutDir/$gene_class/low/0_fold/bootstrap/bootstrap$boot/0_fold.$gene_class.low.omega.txt  >> $OutDir/$gene_class/4_fold_compare/4_fold/bootstrap/bootstrap$boot/4_fold.$gene_class.4_fold_compare.omega.txt

done
##########################################################################################
#####Step11: In order to compare the patterns of purifying and positive selection for different gene classes, in addition to compare their 4-fold synonymous sites patterns, now, I use 4-fold synonymous sites in low-level as the "neutral" ones to compare both high and low gene classes

elif [ "$step" == "11" ]; then

OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out"

for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for anno in {0_fold,utr3,utr5,intron,upstream,downstream}
do
######first creating the relevant folders
Out_anno_real="$OutDir/$gene_class/low_neutral/$anno/real"

if [ ! -d "$Out_anno_real" ]; then
mkdir -p $Out_anno_real
fi
done
done


for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for anno in {0_fold,utr3,utr5,intron,upstream,downstream}
do
for bootstrap_n in {1..200}
do

Out_anno_boot="$OutDir/$gene_class/low_neutral/$anno/bootstrap/bootstrap$bootstrap_n"

if [ ! -d "$Out_anno_boot" ]; then
mkdir -p $Out_anno_boot
fi
done
done
done

#################copy the relevant file to the corresponding folders
###create the .real.txt file
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for anno in {0_fold,utr3,utr5,intron,upstream,downstream}
do
##anno.connectivity.low_neutral.real.txt
sed '4d' $OutDir/$gene_class/high/$anno/real/$anno.$gene_class.high.real.txt > $OutDir/$gene_class/low_neutral/$anno/real/$anno.$gene_class.low_neutral.real.txt
#echo -e "\n" >> $OutDir/$gene_class/low_neutral/$anno/real/$anno.$gene_class.low_neutral.real.txt
tail -n 1 $OutDir/$gene_class/low/$anno/real/$anno.$gene_class.low.real.txt >> $OutDir/$gene_class/low_neutral/$anno/real/$anno.$gene_class.low_neutral.real.txt
#sed '4d' $OutDir/$gene_class/low_neutral/$anno/real/$anno.$gene_class.low_neutral.real.txt > temp && mv temp $OutDir/$gene_class/low_neutral/$anno/real/$anno.$gene_class.low_neutral.real.txt
done
done
####################################################################
###create the .omega.txt file
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for anno in {0_fold,utr3,utr5,intron,upstream,downstream}
do
##$anno.connectivity.$anno_compare.real.txt
sed '2d' $OutDir/$gene_class/high/$anno/real/$anno.$gene_class.high.omega.txt > $OutDir/$gene_class/low_neutral/$anno/real/$anno.$gene_class.low_neutral.omega.txt
tail -n 1 $OutDir/$gene_class/low/$anno/real/$anno.$gene_class.low.omega.txt  >> $OutDir/$gene_class/low_neutral/$anno/real/$anno.$gene_class.low_neutral.omega.txt

done
done

#######################################################################
####for bootstrap dataset
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for anno in {0_fold,utr3,utr5,intron,upstream,downstream}
do
for boot in {1..200}
do
##anno.connectivity.low_neutral.real.txt
sed '4d' $OutDir/$gene_class/high/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.high.real.txt > $OutDir/$gene_class/low_neutral/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.low_neutral.real.txt
tail -n 1 $OutDir/$gene_class/low/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.low.real.txt >> $OutDir/$gene_class/low_neutral/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.low_neutral.real.txt
done
done
done
####################################################################
###create the .omega.txt file
for gene_class in {connectivity,core_gene,eGene,exp_level,exp_variance}   ###I create both low and high even for egene and core_gene
###eGene:high; non-eGene: low
###core_gene: high; non_core: low
do
for anno in {0_fold,utr3,utr5,intron,upstream,downstream}
do
for boot in {1..200}
do
##$anno.connectivity.$anno_compare.bootstrap/bootstrap$boot.txt
sed '2d' $OutDir/$gene_class/high/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.high.omega.txt > $OutDir/$gene_class/low_neutral/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.low_neutral.omega.txt
tail -n 1 $OutDir/$gene_class/low/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.low.omega.txt  >> $OutDir/$gene_class/low_neutral/$anno/bootstrap/bootstrap$boot/$anno.$gene_class.low_neutral.omega.txt

done
done
done








fi



