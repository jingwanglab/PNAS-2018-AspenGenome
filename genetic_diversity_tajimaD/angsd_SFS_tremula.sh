#! /bin/bash -l


#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -n 4
#SBATCH -o angsd_foldedSFS_samtools_24tremula.out
#SBATCH -e angsd_foldedSFS_samtools_24tremula.err
#SBATCH -J angsd_foldedSFS_samtools_24tremula.job
#SBATCH -t 1-00:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL

angsd="/proj/b2011141/tools/angsd0.902/angsd/angsd"
realSFS="/proj/b2011141/tools/angsd0.902/angsd/misc/realSFS"
thetaStat="/proj/b2011141/tools/angsd0.902/angsd/misc/thetaStat"
bam_list_tremula="/proj/b2011141/nobackup/alignments/bwa-mem/asp201/tremula/deduplication/realign/trans/tremula.bam.list"
ref="/proj/b2011141/nobackup/reference/asp201/Potra01-genome.fa"
OutDir="/proj/b2010014/nobackup/population_genetics/tremula/ANGSD/SFS_angsd0.902"
region="/proj/b2010014/nobackup/population_genetics/tremula/rm_cov_MQ0/24tremula.rm.cov.region"

nInd=$(cat $bam_list_tremula | wc -l)
#nChrom=$(echo "2*$nInd" | bc)
nChrom=$nInd

echo $nInd
echo $nChrom


if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

#first generate .saf file
#$angsd -bam $bam_list_tremula -minMapQ 30 -minQ 20 -GL 1 -doSaf 1 -HWE_pval -out $OutDir/tremula_24samples -anc $ref -rf $region -fold 1

#use emOptim2 to optimization
$realSFS $OutDir/tremula_24samples.saf.idx -maxIter 100 -P 4 > $OutDir/tremula_24samples.sfs

#calculate thetas
#$angsd -bam $bam_list_tremula -out $OutDir/tremula_24samples -doThetas 1 -GL 1 -doSaf 1 -anc $ref -rf $chr -pest $OutDir/tremula_24samples.sfs -minMapQ 30 -minQ 20

#calculate Tajimas
#$thetaStat make_bed $OutDir/tremula_24samples.thetas.gz
#$thetaStat do_stat $OutDir/tremula_24samples.thetas.gz -nChr $nChrom
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 10000 -step 10000 -outnames $OutDir/tremula_$1.thetas10kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 50000 -step 50000 -outnames $OutDir/tremula_$1.thetas50kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 100000 -step 100000 -outnames $OutDir/tremula_$1.thetas100kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 500000 -step 500000 -outnames $OutDir/tremula_$1.thetas500kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 1000000 -step 1000000 -outnames $OutDir/tremula_$1.thetas1Mwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 10000 -step 5000 -outnames $OutDir/tremula_$1.thetas10kbwindow.5kbsliding.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 100000 -step 10000 -outnames $OutDir/tremula_$1.thetas100kbwindow.10kbsliding.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 500000 -step 50000 -outnames $OutDir/tremula_$1.thetas500kbwindow.50kbsliding.gz

#calculate posterior probabilities of sample allele frequencies
#$angsd -bam $bam_list_tremula -GL 1 -doSaf 1 -anc $ref -rf $chr -minMapQ 30 -minQ 20 -pest $OutDir/tremula_$1.sfs -out $OutDir/tremula_$1.rf


