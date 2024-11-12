# GOAL: take training fold data, split into 80/20, train Betas in 80, train mixing weights in 20, save linear combination of effect sizes

NUM_THREADS=$(nproc)
TRAINING_PATH=$1 # path + prefix to training data (includes bed/bim/fam and _P.txt and _MG.txt )
ORIGINAL_GENOTYPE_COUNTS=$2 # path to original genotype plink files (where allele counts file is stored)
NUM_PHENOS=$3
thresholds=(0.000005 0.001 0.05)

CV_SUM_dir=${TRAINING_PATH}cv # create directory to store nested cv files
mkdir ${CV_SUM_dir}
CV_SUM_TRAIN=${CV_SUM_dir}/traintrain # prefix for training fold of training fold = used to train Betas 
CV_SUM_TEST=${CV_SUM_dir}/traintest # prefix for testing fold of training fold = used to train mixing weights across PRS

# define number of individuals in the train-train split and train-test split
NUM_TRAINING_FIFTH=$(( $(wc -l < ${TRAINING_PATH}.fam) / 5 ))
head -n $NUM_TRAINING_FIFTH ${TRAINING_PATH}.fam > ${CV_SUM_TEST}.IDS # fifth of the training fold

# initalize files for PRS scores for each phenotype before taking linear combination
for thresh in "${thresholds[@]}"; do
    for pheno in $(seq 1 $NUM_PHENOS); do
        echo -e "FID\tIID\tPRS" > ${CV_SUM_TEST}_prs_predictions_pheno${pheno}_thresh${thresh}_P.txt
    done
done

# Create genotype file for traintrain fold and traintest fold
plink2 \
     --bfile ${TRAINING_PATH} \
     --remove ${CV_SUM_TEST}.IDS \
     --make-bed \
     --out ${CV_SUM_TRAIN} \
     --threads ${NUM_THREADS}
plink2 \
     --bfile ${TRAINING_PATH} \
     --keep ${CV_SUM_TEST}.IDS \
     --make-bed \
     --out ${CV_SUM_TEST} \
     --threads ${NUM_THREADS}

# Filter phenotypes to traintrain fold
python3 ../filter_phenotype.py $TRAINING_PATH $CV_SUM_TRAIN

# GWAS on traintrain fold 
plink2 \
    --bfile $CV_SUM_TRAIN \
    --read-freq ${ORIGINAL_GENOTYPE_COUNTS}'.acount' \
    --pheno ${CV_SUM_TRAIN}'_P.txt' \
    --glm allow-no-covars \
    --out ${CV_SUM_TRAIN}'_P' \
    --threads ${NUM_THREADS}

for pheno in $(seq 1 $NUM_PHENOS); do
    awk 'BEGIN { OFS="\t" } {print $3,$15}' ${CV_SUM_TRAIN}'_P.PHEN'${pheno}'.glm.linear' > ${CV_SUM_TRAIN}_P_SNP.pvalue.PHEN${pheno}
    # plink commands to build PRS score for each phenotype
    plink2 \
        --bfile ${CV_SUM_TEST} \
        --score ${CV_SUM_TRAIN}'_P.PHEN'${pheno}'.glm.linear' 3 4 12 header \
        --q-score-range range_list ${CV_SUM_TRAIN}_P_SNP.pvalue.PHEN${pheno} \
        --out ${CV_SUM_TEST}'_P_PHEN'${pheno} \
        --threads ${NUM_THREADS}
    # Append PRS results for current fold
    for thresh in "${thresholds[@]}"; do
        awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $6}' ${CV_SUM_TEST}'_P'_PHEN${pheno}.${thresh}.sscore >> ${CV_SUM_TEST}_prs_predictions_pheno${pheno}_thresh${thresh}_P.txt
    done
done

# compute a linear combination of PRS to predict maxh phenotype
python3 PRS_sum_fit.py $TRAINING_PATH ${CV_SUM_TEST}_prs_predictions_pheno ${CV_SUM_TEST}.fam ${CV_SUM_TRAIN}'_P' $NUM_PHENOS

