#! /bin/bash -l

#SBATCH -A snic2016-7-89
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o SnIPRE.SwAsp.94samples.input.out
#SBATCH -e SnIPRE.SwAsp.94samples.input.err
#SBATCH -J SnIPRE.SwAsp.94samples.input.job
#SBATCH -t 3-00:00:00


#########The main aim of this script is to create the input file for SnIPRE to run on the 94 SwAsp reseq samples
##Basic steps:  1. Focus on the 22,306 expressed genes from Plos Genetics (2017) paper
#		2. Output file:  gene name, polymorphic 0_fold sites, fixed 0_fold sites, total 0_fold sites, polymorphic 4_fold sites, fixed 4_fold sites, total 4_fold sites were needed
#		3. polymorhpic sites were taken from 94 SwAsp vcf files, and fixed sites were taken from vcf of P.trichocarpa.



module load bioinfo-tools
module load BEDTools
module load R

tri_0_fold="/proj/uppstore2017142/b2010014_nobackup/population_genetics/DFE_expression/outgroup_vcf/trichocarpa/bed_intersect/trichocarpa.frq.filter.0_fold.bed"
tri_4_fold="/proj/uppstore2017142/b2010014_nobackup/population_genetics/DFE_expression/outgroup_vcf/trichocarpa/bed_intersect/trichocarpa.frq.filter.4_fold.bed"
poly_0_fold="/proj/uppstore2017142/b2010014_nobackup/population_genetics/DFE_expression/vcf_count_94samples/anno/SwAsp_94samples.polymorphic.maf.0_fold.bed"
poly_4_fold="/proj/uppstore2017142/b2010014_nobackup/population_genetics/DFE_expression/vcf_count_94samples/anno/SwAsp_94samples.polymorphic.maf.4_fold.bed"
gene_all="/proj/uppstore2017142/b2010014_nobackup/population_genetics/DFE_expression/plos_genetics/S2_file_all_gene_statistics.tsv"
gene_express="/proj/uppstore2017142/b2010014_nobackup/population_genetics/DFE_expression/plos_genetics/S2_file_network_eqtl_gene_statistics.tsv"

###gene, poly_0_fold,fixed_0_fold,total_0_fold,poly_4_fold,fixed_4_fold,total_4_fold

OutDir="/proj/uppstore2017142/b2010014_nobackup/population_genetics/DFE_expression/SnIPRE"

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi


step=$1  ##indicates which step to perfrom
OutDir_poly="$OutDir/tri_poly"  ###the polymorphic sites based on the avaialbe of outgroup sequence

if [ ! -d "$OutDir_poly" ]; then
mkdir -p $OutDir_poly
fi

###Step1. only consider the polymorhic sites that with outgroup species (P.trichocarpa) information available
if [ "$step" == "1" ]; then
bedtools intersect -a $poly_0_fold -b $tri_0_fold -wa > $OutDir_poly/SwAsp_94samples.polymorphic.maf.0_fold.tri.bed 
bedtools intersect -a $poly_4_fold -b $tri_4_fold -wa > $OutDir_poly/SwAsp_94samples.polymorphic.maf.4_fold.tri.bed 

elif [ "$step" == "2" ]; then
###Step2: output the files with columns gene, poly_0_fold,fixed_0_fold,total_0_fold,poly_4_fold,fixed_4_fold,total_4_fold
##for 22,306 expressed genes from Plos Genet paper
poly_0_fold_new=$OutDir_poly/SwAsp_94samples.polymorphic.maf.0_fold.tri.bed
poly_4_fold_new=$OutDir_poly/SwAsp_94samples.polymorphic.maf.4_fold.tri.bed

echo -e "gene\tpoly_0_fold\tfixed_0_fold\ttotal_0_fold\tpoly_4_fold\tfixed_4_fold\ttotal_4_fold" > $OutDir/SwAsp_94samples.SnIPRE.input.txt

for n in {1..22306}
do
gene=$(sed '1d' $gene_express |cut -f 1 |head -n $n |tail -n 1 )
#0_fold
n_poly_0_fold=$(grep "$gene" $poly_0_fold_new |wc -l)
n_total_0_fold=$(grep "$gene" $tri_0_fold |wc -l)
n_fixed_0_fold=$(grep "$gene" $tri_0_fold |awk '$4==2' |awk '$5==0' |wc -l)
#4_fold
n_poly_4_fold=$(grep "$gene" $poly_4_fold_new |wc -l)
n_total_4_fold=$(grep "$gene" $tri_4_fold |wc -l)
n_fixed_4_fold=$(grep "$gene" $tri_4_fold |awk '$4==2' |awk '$5==0' |wc -l)

echo -e "$gene\t$n_poly_0_fold\t$n_fixed_0_fold\t$n_total_0_fold\t$n_poly_4_fold\t$n_fixed_4_fold\t$n_total_4_fold" >> $OutDir/SwAsp_94samples.SnIPRE.input.txt
done


fi





