#! /bin/bash -l

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o est_dfe.expression.bootstrap.SwAsp94.out
#SBATCH -e est_dfe.expression.bootstrap.SwAsp94.err
#SBATCH -J est_dfe.expression.bootstrap.SwAsp94.job
#SBATCH -t 24:00:00



est_dfe="/proj/b2011141/tools/dfe-alpha-release-2.15/est_dfe"
est_alpha_omega="/proj/b2011141/tools/dfe-alpha-release-2.15/est_alpha_omega"
prop_muts_in_s_ranges="/proj/b2011141/tools/dfe-alpha-release-2.15/prop_muts_in_s_ranges"
class=$1 ###connectivity,core_gene,eGene,exp_level,exp_variance
level=$2  ##low or high
anno=$3  ###0_fold,4_fold,utr3,utr5,intron,upstream,downstream
step=$4  ###indicates which step to perform


if [ "$step" == "1" ]; then
#####Step1:site_class_0
###submit jobs
##for file in {connectivity,core_gene,eGene,exp_level,exp_variance}; do for file2 in {low,high}; do for file3 in {0_fold,4_fold,intron,utr3,utr5,upstream,downstream}; do sbatch est_dfe.expression.bootstrap.SwAsp94.sh $file $file2 $file3 1; done; done; done
for bootstrap in {1..200}
do

OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/$class/$level/$anno/bootstrap/bootstrap$bootstrap"

echo  "data_path_1    /proj/b2011141/tools/dfe-alpha-release-2.15/data" > $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt 
echo  "data_path_2    /proj/b2011141/tools/dfe-alpha-release-2.15/data/data-three-epoch" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "sfs_input_file  $OutDir/$anno.$class.$level.real.txt" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "est_dfe_results_dir   $OutDir/results_dir_neut/" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "site_class  0" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "fold   1">> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "epochs  2" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "search_n2  1" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "t2_variable   1" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt
echo "t2 50" >> $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt

$est_dfe -c $OutDir/est_dfe.$class.$level.$anno.site_class-0.txt

done

#####Step2:site_class_1
elif [ "$step" == "2" ]; then
###submit jobs
##for file in {connectivity,core_gene,eGene,exp_level,exp_variance}; do for file2 in {low,high}; do for file3 in {0_fold,4_fold,intron,utr3,utr5,upstream,downstream}; do sbatch est_dfe.expression.SwAsp94.sh $file $file2 $file3 2; done; done; done

for bootstrap in {1..200}
do

OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/$class/$level/$anno/bootstrap/bootstrap$bootstrap"

echo  "data_path_1    /proj/b2011141/tools/dfe-alpha-release-2.15/data" > $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo  "data_path_2    /proj/b2011141/tools/dfe-alpha-release-2.15/data/data-three-epoch" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "sfs_input_file  $OutDir/$anno.$class.$level.real.txt" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "est_dfe_results_dir   $OutDir/results_dir_sel/" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "est_dfe_demography_results_file   $OutDir/results_dir_neut/est_dfe.out" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "site_class  1" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "fold   1">> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "epochs  2" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "mean_s_variable  1" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "mean_s   -0.1" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "beta_variable  1" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt
echo "beta   0.5" >> $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt

$est_dfe -c $OutDir/est_dfe.$class.$level.$anno.site_class-1.txt

$prop_muts_in_s_ranges -c $OutDir/results_dir_sel/est_dfe.out -o $OutDir/results_dir_sel/prop_muls_in_s_ranges.output

done

#####Step3: est_alpha_omega
elif [ "$step" == "3" ]; then
###submit jobs
##for file in {connectivity,core_gene,eGene,exp_level,exp_variance}; do for file2 in {low,high}; do for file3 in {0_fold,4_fold,intron,utr3,utr5,upstream,downstream}; do sbatch est_dfe.expression.SwAsp94.sh $file $file2 $file3 3; done; done; done

for bootstrap in {1..200}
do
OutDir="/proj/b2010014/nobackup/population_genetics/DFE_expression/Out/$class/$level/$anno/bootstrap/bootstrap$bootstrap"

echo "data_path_1    /proj/b2011141/tools/dfe-alpha-release-2.15/data/" > $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt
echo "divergence_file  $OutDir/$anno.$class.$level.omega.txt" >> $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt
echo "est_alpha_omega_results_file   $OutDir/est_alpha_omega.$anno.$class.$level.out" >> $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt
echo "est_dfe_results_file  $OutDir/results_dir_sel/est_dfe.out" >> $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt
echo "neut_egf_file  $OutDir/results_dir_neut/neut_egf.out" >> $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt
echo "sel_egf_file  $OutDir/results_dir_sel/sel_egf.out" >> $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt
echo "do_jukes_cantor  1" >> $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt
echo "remove_poly  0" >> $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt

$est_alpha_omega -c $OutDir/est_dfe.$class.$level.$anno.alpha_omega.txt

done
fi



