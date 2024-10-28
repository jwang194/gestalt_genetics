#!/bin/bash
#$ -o /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/logs/gwas_phq.o$TASK_ID
#$ -e /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/logs/gwas_phq.e$TASK_ID
#$ -l h_data=8G
#$ -l time=47:59:00
#$ -l highp
#$ -pe shared 4
#$ -t 1-17
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

. /u/local/Modules/default/init/modules.sh
module load zlib

PN=$(sed "${SGE_TASK_ID}q;d" /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/data/phens/phq_colnames.txt)

for CHR in {1..22}
do
    /u/home/r/ribo7412/project-zaitlenlab/bin/plink2_avx2 \
        --bfile  /u/project/sgss/UKBB/data/geno_QC/impSNPs_unrel_EUR_maf01_info.9/impSNPs_unrel_EUR_maf01_info.9_geno.1_hwe1em7_chr"$CHR" \
        --read-freq  /u/project/sgss/UKBB/data/geno_QC/impSNPs_unrel_EUR_maf01_info.9/impSNPs_unrel_EUR_maf01_info.9_geno.1_hwe1em7_chr"$CHR".afreq \
        --covar /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/data/covs/PC1_5_Age_Sex.qcovar \
        --pheno /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/data/phens/phq.pheno \
        --pheno-name ${PN} \
        --threads 4 \
        --memory 127000 \
        --vif 9999 \
        --covar-variance-standardize \
        --glm hide-covar \
        --out /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/results/sumstats/phq_${CHR}
done
