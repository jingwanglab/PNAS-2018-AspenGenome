#! /bin/bash -l

#SBATCH -A b2011141
#SBATCH -p core
#SBATCH -o plink_make_bed.out
#SBATCH -e plink_make_bed.err
#SBATCH -J plink_make_bed.job
#SBATCH -t 1-00:00:00
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

make_bed=$OutDir/make_bed

if [ ! -d "$make_bed" ]; then
mkdir -p $make_bed
fi

$plink --ped $ped --map $map --hwe 0.0001 --make-bed --out $make_bed/$Output --recode
