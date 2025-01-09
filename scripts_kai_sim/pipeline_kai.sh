#!/bin/bash
#$ -o /u/scratch/k/kaia/JobLogs/gestalt/pipeline.o$JOB_ID_$TASK_ID
#$ -e /u/scratch/k/kaia/JobLogs/gestalt/pipeline.e$JOB_ID_$TASK_ID
#$ -l h_data=128G
#$ -l time=2:00:00
#$ -t 1-20
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

R=$SGE_TASK_ID
#R=1
N=${1}
M=${2} # total number of snps (includes both causal and noncausal)
B=${3}
C=${4}
P_DIR=${5} # directory to save phenotype files (Jeremy = "../data/sim/phenos/", Kai = "../../data/sim/phenos/")
GWAS_DIR=${6} # directory to save GWAS output (Jeremy = "../gwas/sim/", Kai = "../../data/sim/gwas/")
SCRATCH_DIR=${7}

. /u/local/Modules/default/init/modules.sh
#. /u/home/j/jwang194/.profile

module load mamba
mamba activate gestalt

#pyl torch

#cd /u/home/j/jwang194/zp/zaitlen/gestalt/scripts
cd /u/home/k/kaia/GESTALT/gestalt_genetics/scripts_kai_sim

python3 src/generate.py $N $M ${SCRATCH_DIR}'configs/'${C}'/'${C}'_assignments.txt' ${SCRATCH_DIR}'configs/'${C}'/'${C}'_gg.txt' ${SCRATCH_DIR}'configs/'${C}'/'${C}'_ge.txt' $B $R $C $P_DIR

#condaload maxgcp

python3 src/mgp.py ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_P.txt' ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_PE.txt' ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_MG.txt' ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_SD.txt'

python3 src/rename.py ${N} ${M} ${C} ${R} ${P_DIR}

source ~/.bashrc

plink2 \
    --bfile $B \
    --read-freq $B'.acount' \
    --pheno ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_P.txt' \
    --glm allow-no-covars \
    --out ${GWAS_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_P'

#--snps "SNP_0-SNP_"$(($M-1)) \
    # before we set snps 0-M as causal 
    # now snps are chosen randomly, and the true causal variants can be found in configs variant assignment file 
    # bfiles filtered for these causal snps should be input to $B

plink2 \
    --bfile $B \
    --read-freq $B'.acount' \
    --pheno ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_MG.txt' \
    --glm allow-no-covars \
    --out ${GWAS_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_MG'

# --snps "SNP_0-SNP_"$(($M-1)) \

plink2 \
    --bfile $B \
    --read-freq $B'.acount' \
    --pheno ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_SD.txt' \
    --glm allow-no-covars \
    --out ${GWAS_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_SD'

# --snps "SNP_0-SNP_"$(($M-1)) \


