#!/bin/bash
#$ -o /u/home/k/kaia/JobLogs/gestalt/pipeline.o$JOB_ID_$TASK_ID
#$ -e /u/home/k/kaia/JobLogs/gestalt/pipeline.e$JOB_ID_$TASK_ID
#$ -l h_data=64G
#$ -l time=16:00:00
#$ -t 1-25
# Email address to notify
#$ -M kakamatsu@g.ucla.edu
# Notify when
#$ -m bea
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

R=$SGE_TASK_ID
N=${1}
M=${2}
B=${3}
C=${4}
P_DIR=${5} # directory to save phenotype files (Jeremy = "../data/sim/phenos/", Kai = "../../data/sim/phenos/")
GWAS_DIR=${6} # directory to save GWAS output (Jeremy = "../gwas/sim/", Kai = "../../data/sim/gwas/")

. /u/local/Modules/default/init/modules.sh
#. /u/home/j/jwang194/.profile

module load mamba
mamba activate gestalt

#pyl torch

#cd /u/home/j/jwang194/zp/zaitlen/gestalt/scripts
cd /u/home/k/kaia/GESTALT/gestalt_genetics/scripts_kai

python3 generate.py $N $M 'configs/'${C}'_assignments.txt' 'configs/'${C}'_gg.txt' 'configs/'${C}'_ge.txt' $B $R $C $P_DIR

#condaload maxgcp

python3 mgp.py ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_P.txt' ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_PE.txt' ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_MG.txt' ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_SD.txt'

python3 rename.py ${N} ${M} ${C} ${R} ${P_DIR}

source ~/.bashrc

plink2 \
    --bfile $B \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq $B'.acount' \
    --pheno ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_P.txt' \
    --glm allow-no-covars \
    --out ${GWAS_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_P'

plink2 \
    --bfile $B \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq $B'.acount' \
    --pheno ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_MG.txt' \
    --glm allow-no-covars \
    --out ${GWAS_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_MG'

plink2 \
    --bfile $B \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq $B'.acount' \
    --pheno ${P_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_SD.txt' \
    --glm allow-no-covars \
    --out ${GWAS_DIR}${N}'_'${M}'_'${C}'_rep'${R}'_SD'
