### 1. loop through all configurations ###
### 2. submit jobs to slurm ###

# phenotypes should already be simulated 
    # see: /u/home/k/kaia/GESTALT/gestalt_genetics/scripts_kai_sim

. /u/local/Modules/default/init/modules.sh
module load mamba
mamba activate gestalt
cd /u/home/k/kaia/GESTALT/gestalt_genetics/PRS/simulations

PLINK2_PATH=/u/home/k/kaia/STATSGEN/plink/plink2
sim_dir=/u/home/${USER:0:1}/${USER}/GESTALT/sim
scratch_sim_dir=/u/scratch/${USER:0:1}/${USER}/GESTALT/sim 
    # all files saved in scratch to conserve disk space
    # make sure to save gwas power results and PRS results to home so they are not deleted

polygenicity=(0.01 0.05 0.1 0.2)
shared_fraction=(0.05 0.2 0.4)
genetic_covariance=(0.01 0.04 0.08) # of shared components (total genetic covariance between traits is 1/2 of this value)
environmental_covariance=(-0.1 0.1)

### incomplete CV run, switch to predicting on held-out cohort because CV takes too long ### 
#polygenicity=(0.01)
#shared_fraction=(0.4)
#genetic_covariance=(0.01 0.04 0.08) # of shared components (total genetic covariance between traits is 1/2 of this value)
#environmental_covariance=(-0.1 0.1)

#polygenicity=(0.05 0.1 0.2)
#shared_fraction=(0.05 0.2 0.4)
#genetic_covariance=(0.01 0.04 0.08) # of shared components (total genetic covariance between traits is 1/2 of this value)
#environmental_covariance=(-0.1 0.1)
### ---

num_traits=10
N=50000
M=98163

#genotypes=${sim_dir}/genotypes/N50kM1m/sim_N50kM1m # prefix to plink genotype files (full, not just causal) 
    # gwas with 1 million and 500k variants cannot run with 128GB of memory
#genotypes=${sim_dir}/genotypes/N50kM100k_LD/sim_N50kM100k 
#genotypes=${sim_dir}/genotypes/N1kM10k/sim_N1000_M10000
genotypes=${sim_dir}/genotypes/N50kM100k_LD/maf_0.01_filter/sim_N50kM100k_maf

## to improve runtime, change LD reference to a subset of the population ##
    # --clump on a large dataset can take a long time
ldref=${sim_dir}/genotypes/N50kM100k_LD/maf_0.01_filter/ldref_subset/sim_N50kM100k_maf_subset
#ldref=${sim_dir}/genotypes/N1kM10k/sim_N1000_M10000

for poly in "${polygenicity[@]}"; do
    for shared in "${shared_fraction[@]}"; do
        for cov in "${genetic_covariance[@]}"; do
            for env in "${environmental_covariance[@]}"; do
                echo "Polygenicity: $poly, Shared causal variants: $shared, Genetic covariance: $cov, Environmental covariance: $env"
                setting=traits_${num_traits}_causal_${poly}_shared_${shared}_uniform_rg_${cov}_random_re_${env}
                phenotype_dir=${scratch_sim_dir}/phenos/${setting}
                out_dir=${scratch_sim_dir}/PRS/${setting}
                mkdir -p $out_dir
                echo "Running pipeline"
                qsub prs_pipeline.sh $N $M $setting ${genotypes} ${ldref} ${phenotype_dir} ${out_dir} ${PLINK2_PATH}
            done
        done
    done
done