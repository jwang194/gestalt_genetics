#!/bin/bash
#$ -o /u/scratch/k/kaia/JobLogs/gestalt/PRSsum.o$JOB_ID
#$ -e /u/scratch/k/kaia/JobLogs/gestalt/PRSsum.e$JOB_ID
#$ -l h_data=64G
#$ -l time=2:00:00
#$ -l arch="amd-epyc-7642|i-E5-2670v3|intel-gold*"

. /u/local/Modules/default/init/modules.sh
module load mamba
mamba activate gestalt
cd /u/home/k/kaia/GESTALT/gestalt_genetics/scripts_kai/PRS/PRS_updated_pipeline/

# iterate through simulation replicates
for i in {1..25}; do
    bash PRS_plinkCV.sh $i
done
#for i in {1..25}; do
#    bash PRS_plinkCV2.sh $i
#done
#for i in {1..25}; do
#    bash PRS_plinkCV3.sh $i
#done
#for i in {1..25}; do
#    bash PRS_plinkCV4.sh $i
#done

