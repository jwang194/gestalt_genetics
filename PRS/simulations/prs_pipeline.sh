#!/bin/bash
#$ -o /u/scratch/k/kaia/JobLogs/gestalt/PRS.o$JOB_ID_$TASK_ID
#$ -e /u/scratch/k/kaia/JobLogs/gestalt/PRS.e$JOB_ID_$TASK_ID
#$ -l h_data=128G
#$ -l time=3:00:00
#$ -t 1-20
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

R=$SGE_TASK_ID
#R=1
N=${1} # number of individuals
M=${2} # total number of snps (includes both causal and noncausal)
SETTING=${3}
GENOTYPES=${4} # path + prefix to plink genotype files
LDREF=${5}
P_DIR=${6} 
OUT_DIR=${7} 
PLINK2_PATH=${8} # path to plink2 executable

. /u/local/Modules/default/init/modules.sh
module load mamba
mamba activate gestalt
cd /u/home/k/kaia/GESTALT/gestalt_genetics/PRS/simulations


echo 'PRS for replicate' $R
#CV_DIR=${OUT_DIR}/cv_rep${R}
#mkdir -p $CV_DIR
INTERMEDIATE_DIR=${OUT_DIR}/intermediate_rep${R}
mkdir -p $INTERMEDIATE_DIR
P_FILE=${P_DIR}/${N}_${M}_${SETTING}_rep${R} # prefix to _P.txt and _MG.txt phenotype files
OUT_FILE=${OUT_DIR}/${N}_${M}_${SETTING}_rep${R}_prs # prefix to output files
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
    # number of threads
    # evaluate R-squared between true and predicted phenotypes
### P+T details ###
    # prune using --clump-p1 1 --clump-r2 0.1 --clump-kb 250 on plink
    # threshold at various p-values, evaluate on validation cohort, and choose the best one [5e-8 5e-6 5e-4 0.05]
python3 src/PRS.py --plink2_path $PLINK2_PATH \
                --genotypes $GENOTYPES \
                --ldref $LDREF \
                --phenotype_prefix $P_FILE \
                --intermediate_dir $INTERMEDIATE_DIR \
                --output_prefix $OUT_FILE \
                --model_type P,MG,SUMPRS_MG,PRSSUM,SUMPRS \
                --threads 8 \
                --evaluate_r2 \
                --p_values 5e-8,5e-6,5e-4,0.05 \
                --validation_test_fraction 0.1

#python3 src/PRS_CV.py --plink2_path $PLINK2_PATH \
#                --genotypes $GENOTYPES \
#                --ldref $LDREF \
#                --phenotype_prefix $P_FILE \
#                --intermediate_dir $CV_DIR \
#                --output_prefix $OUT_FILE \
#                --model_type P,MG,SUMPRS_MG,PRSSUM,SUMPRS \
#                --cv_folds 5 \
#                --threads 8 \
#                --evaluate_r2 \
#                --p_values 5e-8,5e-6,5e-4,0.05 \
#                --validation_fraction 0.1
