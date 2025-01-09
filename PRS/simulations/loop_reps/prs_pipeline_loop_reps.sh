#!/bin/bash
#$ -o /u/scratch/k/kaia/JobLogs/gestalt/PRS.o$JOB_ID
#$ -e /u/scratch/k/kaia/JobLogs/gestalt/PRS.e$JOB_ID
#$ -l h_data=128G
#$ -l time=20:00:00
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

N=${1} # number of individuals
M=${2} # total number of snps (includes both causal and noncausal)
declare -i REPS=${3} # total number of replicates 
SETTING=${4}
GENOTYPES=${5} # path + prefix to plink genotype files
LDREF=${6}
P_DIR=${7} 
OUT_DIR=${8} 
PLINK2_PATH=${9} # path to plink2 executable

. /u/local/Modules/default/init/modules.sh
module load mamba
mamba activate gestalt
cd /u/home/k/kaia/GESTALT/gestalt_genetics/PRS/simulations

# for a given config, iteration through all replicates
for i in $(seq 1 $REPS); do
    echo 'PRS for replicate' $i
    CV_DIR=${OUT_DIR}/cv_rep${i}
    mkdir -p $CV_DIR
    P_FILE=${P_DIR}/${N}_${M}_${SETTING}_rep${i} # prefix to _P.txt and _MG.txt phenotype files
    OUT_FILE=${OUT_DIR}/${N}_${M}_${SETTING}_rep${i}_prs # prefix to output files
    # python script to evaluate PRS, traditional pruning and thresholding method
    ### INPUTS ###
        # path to plink2 executable
        # plink genotype files 
        # phenotype file prefix
        # intermediate directory for cross-validation files 
        # output file prefix
    ### OPTIONS ### 
        # model/pheno type: comma-separated list from ['P', 'MG', 'SUMPRS_MG', 'PRSSUM', 'SUMPRS']
            # P, PRSSUM, SUMPRS requires _P.txt file
            # MG, SUMPRS_MG requires _MG.txt file
        # number of CV folds 
        # number of threads
        # evaluate R-squared between true and predicted phenotypes
    ### P+T details ###
        # prune using --clump-p1 1 --clump-r2 0.1 --clump-kb 250 on plink
        # threshold at various p-values, evaluate via nested cross validation, and choose the best one [5e-8 5e-6 5e-4 0.05]
    python3 src/PRS.py --plink2_path $PLINK2_PATH \
                    --genotypes $GENOTYPES \
                    --ldref $LDREF \
                    --phenotype_prefix $P_FILE \
                    --intermediate_dir $CV_DIR \
                    --output_prefix $OUT_FILE \
                    --model_type P,MG,SUMPRS_MG,PRSSUM,SUMPRS \
                    --cv_folds 5 \
                    --threads 8 \
                    --evaluate_r2 \
                    --p_values 5e-8,5e-6,5e-4,0.05 \
                    --validation_fraction 0.1
done
