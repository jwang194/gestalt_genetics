### 1. define all parameters to test in simulations###
### 2. loop through all possible combinations of parameters, create config files ###
### 3. submit jobs to slurm ###

. /u/local/Modules/default/init/modules.sh
module load mamba
mamba activate gestalt
cd /u/home/k/kaia/GESTALT/gestalt_genetics/scripts_kai_sim

sim_dir=/u/home/${USER:0:1}/${USER}/GESTALT/sim
scratch_sim_dir=/u/scratch/${USER:0:1}/${USER}/GESTALT/sim 
    # all files saved in scratch to conserve disk space
    # make sure to save gwas power results and PRS results to home so they are not deleted

polygenicity=(0.01 0.05 0.1 0.2)
shared_fraction=(0.05 0.2 0.4)
#genetic_covariance=(0.01 0.04 0.08) # of shared components (total genetic covariance between traits is 1/2 of this value)
genetic_covariance=(0.001)
environmental_covariance=(-0.1 0.1)

num_traits=10
N=50000
M=98163


for poly in "${polygenicity[@]}"; do
    for shared in "${shared_fraction[@]}"; do
        for cov in "${genetic_covariance[@]}"; do
            for env in "${environmental_covariance[@]}"; do
                echo "Polygenicity: $poly, Shared causal variants: $shared, Genetic covariance: $cov, Environmental covariance: $env"
                setting=traits_${num_traits}_causal_${poly}_shared_${shared}_uniform_rg_${cov}_random_re_${env}
                mkdir -p ${scratch_sim_dir}/phenos/${setting}/
                mkdir -p ${scratch_sim_dir}/gwas/${setting}/
                # cov is the genetic covariance of the shared components (not the total genetic covariance between traits)
                    # by design, the additive genetic variance of the specific and shared components are both set to the same value, and the betas are scaled so that the trait heritability is what we want
                    # therefore, the total genetic covariance between traits is 1/2 of the genetic covariance of the shared components
                python3 configs/create_configs.py 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 ${sim_dir}/genotypes/N50kM100k_LD/maf_0.01_filter/sim_N50kM100k_maf.bim $cov $env $M $poly $shared ${scratch_sim_dir}/configs # configs files written in scratch
                echo "Finish configs, filtering genotypes to causal variants"
                # we want to simulate phenotypes and test GWAS power so we can filter out all noncausal variants 
                mkdir -p ${scratch_sim_dir}/genotypes/N50kM100k_LD/${setting} # genotypes for each setting in scratch
                plink2 --bfile ${sim_dir}/genotypes/N50kM100k_LD/maf_0.01_filter/sim_N50kM100k_maf --extract ${scratch_sim_dir}/configs/${setting}/${setting}_assignments.txt --make-bed --out ${scratch_sim_dir}/genotypes/N50kM100k_LD/${setting}/${setting}
                plink2 --bfile ${scratch_sim_dir}/genotypes/N50kM100k_LD/${setting}/${setting} --freq counts --out ${scratch_sim_dir}/genotypes/N50kM100k_LD/${setting}/${setting}
                echo "Running pipeline"
                qsub pipeline_kai.sh $N $M ${scratch_sim_dir}/genotypes/N50kM100k_LD/${setting}/${setting} $setting ${scratch_sim_dir}/phenos/${setting}/ ${scratch_sim_dir}/gwas/${setting}/ ${scratch_sim_dir}/
            done
        done
    done
done


