#!/bin/bash
#$ -o /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/logs/gwas_dn_specific.o$TASK_ID
#$ -e /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/logs/gwas_dn_specific.e$TASK_ID
#$ -l h_data=8G
#$ -l time=11:59:00
#$ -l highp
#$ -pe shared 4
#$ -t 3-22
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

. /u/local/Modules/default/init/modules.sh
module load zlib

CHR=${SGE_TASK_ID}
PN=${1}

/u/home/r/ribo7412/project-zaitlenlab/bin/plink2_avx2 \
    --bfile  /u/project/sgss/UKBB/data/geno_QC/impSNPs_unrel_EUR_maf01_info.9/impSNPs_unrel_EUR_maf01_info.9_geno.1_hwe1em7_chr"$CHR" \
    --read-freq  /u/project/sgss/UKBB/data/geno_QC/impSNPs_unrel_EUR_maf01_info.9/impSNPs_unrel_EUR_maf01_info.9_geno.1_hwe1em7_chr"$CHR".afreq \
    --covar /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/data/covs/PC1_5_Age_Sex.qcovar \
    --pheno /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/data/phens/depression_neuroticism.pheno \
    --pheno-name ${PN} \
    --threads 4 \
    --memory 127000 \
    --vif 9999 \
    --covar-variance-standardize \
    --glm hide-covar \
    --out /u/home/j/jwang194/zp/zaitlen/gestalt/gwas/results/sumstats/depr_neur_${CHR}
