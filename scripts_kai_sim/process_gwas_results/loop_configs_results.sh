#!/bin/bash
#$ -o /u/scratch/k/kaia/JobLogs/gestalt/GWAS_results.o$JOB_ID
#$ -e /u/scratch/k/kaia/JobLogs/gestalt/GWAS_results.e$JOB_ID
#$ -l h_data=64G
#$ -l time=3:00:00

. /u/local/Modules/default/init/modules.sh
module load mamba
mamba activate gestalt
cd /u/home/k/kaia/GESTALT/gestalt_genetics/scripts_kai_sim/process_gwas_results

polygenicity=(0.01 0.05 0.1 0.2)
shared_fraction=(0.05 0.2 0.4)
genetic_covariance=(0.001 0.01 0.04 0.08) # of shared components (total genetic covariance between traits is 1/2 of this value)
environmental_covariance=(-0.1 0.1)

num_traits=10
N=50000
M=98163
R=20

BETAS=/u/scratch/k/kaia/GESTALT/sim/phenos
GWAS=/u/scratch/k/kaia/GESTALT/sim/gwas
OUT=/u/home/k/kaia/GESTALT/sim/gwas_results
CONFIGS=/u/scratch/k/kaia/GESTALT/sim/configs

for poly in "${polygenicity[@]}"; do
    for shared in "${shared_fraction[@]}"; do
        for cov in "${genetic_covariance[@]}"; do
            for env in "${environmental_covariance[@]}"; do
                echo "Polygenicity: $poly, Shared causal variants: $shared, Genetic covariance: $cov, Environmental covariance: $env"
                setting=traits_${num_traits}_causal_${poly}_shared_${shared}_uniform_rg_${cov}_random_re_${env}
                BETAS_DIR=${BETAS}/${setting}
                GWAS_DIR=${GWAS}/${setting}
                OUT_DIR=${OUT}/${setting}
                mkdir -p $OUT_DIR
                python3 evaluate.py $BETAS_DIR $GWAS_DIR $OUT_DIR $setting $N $M $num_traits $R $CONFIGS
            done
        done
    done
done