#!/bin/bash
#$ -o /u/home/j/jwang194/zp/zaitlen/gestalt/scripts/logs/pipeline.o$TASK_ID
#$ -e /u/home/j/jwang194/zp/zaitlen/gestalt/scripts/logs/pipeline.e$TASK_ID
#$ -l h_data=127G
#$ -l time=12:00:00
#$ -t 1-25
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

R=$SGE_TASK_ID
N=${1}
M=${2}
B=${3}
C=${4}

. /u/local/Modules/default/init/modules.sh
. /u/home/j/jwang194/.profile

pyl torch

cd /u/home/j/jwang194/zp/zaitlen/gestalt/scripts

python3 generate.py $N $M 'configs/'${C}'_assignments.txt' 'configs/'${C}'_gg.txt' 'configs/'${C}'_ge.txt' $B $R $C

condaload maxgcp

python3 mgp.py '../data/sim/phenos/'${N}'_'${M}'_'${C}'_rep'${R}'_P.txt' '../data/sim/phenos/'${N}'_'${M}'_'${C}'_rep'${R}'_PE.txt' '../data/sim/phenos/'${N}'_'${M}'_'${C}'_rep'${R}'_MG.txt' '../data/sim/phenos/'${N}'_'${M}'_'${C}'_rep'${R}'_SD.txt'

python3 rename.py ${N} ${M} ${C} ${R}

plink2_avx2 \
    --bfile $B \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq $B'.acount' \
    --pheno '../data/sim/phenos/'${N}'_'${M}'_'${C}'_rep'${R}'_P.txt' \
    --glm allow-no-covars \
    --out '../gwas/sim/'${N}'_'${M}'_'${C}'_rep'${R}'_P'

plink2_avx2 \
    --bfile $B \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq $B'.acount' \
    --pheno '../data/sim/phenos/'${N}'_'${M}'_'${C}'_rep'${R}'_MG.txt' \
    --glm allow-no-covars \
    --out '../gwas/sim/'${N}'_'${M}'_'${C}'_rep'${R}'_MG'

plink2_avx2 \
    --bfile $B \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq $B'.acount' \
    --pheno '../data/sim/phenos/'${N}'_'${M}'_'${C}'_rep'${R}'_SD.txt' \
    --glm allow-no-covars \
    --out '../gwas/sim/'${N}'_'${M}'_'${C}'_rep'${R}'_SD'
