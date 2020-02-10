#! /bin/bash -l

#SBATCH -A b2011141
#SBATCH -p core
#SBATCH -o plink_thin_snp.out
#SBATCH -e plink_thin_snp.err
#SBATCH -J plink_thin_snp.job
#SBATCH -t 5:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL

#sbatch plink_thin_snp.sh /proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/trichocarpa/plink/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.trichocarpa.plink.map /proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/trichocarpa/plink/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.trichocarpa.plink.ped

plink="/proj/b2011141/tools/plink_1.9/plink"
map=$1
ped=$2

OutDir=`dirname $1`
Out=${1##*/}
Output=${Out%.map}

###Step 1, because the output map file from vcftools have the Chr information with 0, so we need to first replace the 0 by the real chr number
cut -f 1 -d ":" $map | cut -f 2 | sed 's/Chr0//g' | sed 's/Chr//g' > $OutDir/chr
cut -f 2,3,4 $map > $OutDir/map
paste $OutDir/chr $OutDir/map > $map && rm $OutDir/chr $OutDir/map

thin_ld=$OutDir/thin_ld

if [ ! -d "$thin_ld" ]; then
mkdir -p $thin_ld
fi
##Step 2 filtering based on maf, missing (already done), and hwe then thin the SNPs with only 100kb SNPs left for further analysis

$plink --ped $ped --map $map --hwe 0.0001 --thin-count 100000 --out $thin_ld/$Output --recode
#$plink --ped $ped --map $map --hwe 0.0001 --out $thin_ld/$Output --recode
$plink --ped $thin_ld/$Output.ped --map $thin_ld/$Output.map --r2 --allow-no-sex --ld-window 1000 --ld-window-kb 50 --ld-window-r2 0 --out $thin_ld/$Output.thin.r2

rm $thin_ld/$Output.ped $thin_ld/$Output.map $thin_ld/$Output.nosex $thin_ld/$Output.log
