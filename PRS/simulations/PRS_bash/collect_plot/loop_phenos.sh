PRS_file_dir='/u/home/k/kaia/GESTALT/data/sim/PRS/N50kM10k_10traits_rg0.005_re-0.1_all/' # directory to PRS predictions ex. (/u/home/k/kaia/GESTALT/data/sim/PRS/N50kM10k_5traits_rg0.1_re-0.1/)
Pheno_file_prefix='/u/home/k/kaia/GESTALT/data/sim/phenos/N50kM10k_10traits_rg0.005_re-0.1/50000_10000_traits_10_shared_0.2_specific_0.08_uniform_rg_0.005_random_re_-0.1_rep' # phenotype file dir + prefix ex. (/u/home/k/kaia/GESTALT/data/sim/phenos/N50kM10k_5traits_rg0.1_re-0.1/50000_10000_all_and_ind_overlaps_uniform_gg_random_ge_0.1_rep)
Num_pheno=10
Pheno_types=('P' 'MG' 'MG_sumPRS' 'SUM_sumPRS' 'PRS_SUM') # P or MG or MG_sumPRS or SUM_sumPRS or PRS_SUM
out_dir=/u/home/k/kaia/GESTALT/data/sim/PRS/N50kM10k_10traits_rg0.005_re-0.1_all/results/

for pheno_type in "${Pheno_types[@]}";do 
    Rscript collect_results.R $PRS_file_dir $Pheno_file_prefix $Num_pheno $pheno_type $out_dir
done
