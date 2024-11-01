source ~/.bashrc

rep=$1
GENOTYPE_PREFIX="/u/home/k/kaia/GESTALT/data/sim/genotypes/N50kM100k/sim1" # Path to PLINK binary files without extension
PHENOTYPE_FILE="/u/home/k/kaia/GESTALT/data/sim/phenos/N50kM10k_5traits_rg0.1_re-0.1/" # Path to true phenotype file (1000_10000_testing_rep1234_P.txt, 1000_10000_testing_rep1234_MG.txt)
NUM_FOLDS=5 # Number of folds for cross-validation
OUT_DIR="/u/home/k/kaia/GESTALT/data/sim/PRS/N50kM10k_5traits_rg0.1_re-0.1/" # Output directory
setting=50000_10000_all_and_ind_overlaps_uniform_gg_random_ge_0.1_rep${rep}
M=10000 # number of snps
NUM_PHENOS=5


# Split the fam file into N random parts (testing fold for every iteration)
mkdir ${OUT_DIR}
mkdir ${OUT_DIR}cv
shuf "${GENOTYPE_PREFIX}.fam" | split -l $(( $(wc -l < "${GENOTYPE_PREFIX}.fam") / $NUM_FOLDS )) - "${OUT_DIR}cv/fam_fold_"
num=1
for file in "${OUT_DIR}cv/fam_fold_"*; do
    mv "$file" "${OUT_DIR}cv/fam_fold${num}.IDS" # just renaming here
    ((num++))
done

# P-value thresholds for PRS scores (NOTE: in simulations all snps are independent so no LD pruning)
thresholds=(0.000005 0.001 0.05)
> range_list
for thresh in "${thresholds[@]}"; do
    echo "$thresh 0 $thresh" >> range_list # input to plink later
done

# Initialize file to store final PRS predictions
for pheno in $(seq 1 $NUM_PHENOS); do
    for thresh in "${thresholds[@]}"; do
        echo -e "FID\tIID\tPRS" > ${OUT_DIR}prs_predictions_pheno${pheno}_thresh${thresh}rep${rep}_P.txt
        echo -e "FID\tIID\tPRS" > ${OUT_DIR}prs_predictions_pheno${pheno}_thresh${thresh}rep${rep}_MG.txt
    done
done

# Loop through each fold as test set
for test_fold in $(seq 1 $NUM_FOLDS); do
    echo "Processing fold $test_fold..."

    GENOTYPE_TRAIN=${OUT_DIR}cv/train_fold${test_fold}
    GENOTYPE_TEST=${OUT_DIR}cv/test_fold${test_fold}

    # Create plink files for train and test
    plink2 --bfile $GENOTYPE_PREFIX --remove ${OUT_DIR}cv/fam_fold${test_fold}.IDS --make-bed --out ${GENOTYPE_TRAIN}
    plink2 --bfile $GENOTYPE_PREFIX --keep ${OUT_DIR}cv/fam_fold${test_fold}.IDS --make-bed --out ${GENOTYPE_TEST}

    # Filter phenotype file to keep only training individuals
    python3 filter_phenotype.py ${PHENOTYPE_FILE}${setting} $GENOTYPE_TRAIN

    # --- GWAS on training individuals
    plink2 \
    --bfile $GENOTYPE_TRAIN \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq ${GENOTYPE_PREFIX}'.acount' \
    --pheno ${GENOTYPE_TRAIN}'_P.txt' \
    --glm allow-no-covars \
    --out ${GENOTYPE_TRAIN}'_P'

    plink2 \
    --bfile $GENOTYPE_TRAIN \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --read-freq ${GENOTYPE_PREFIX}'.acount' \
    --pheno ${GENOTYPE_TRAIN}'_MG.txt' \
    --glm allow-no-covars \
    --out ${GENOTYPE_TRAIN}'_MG'
    # ---

    for pheno in $(seq 1 $NUM_PHENOS); do
        awk 'BEGIN { OFS="\t" } {print $3,$15}' ${GENOTYPE_TRAIN}'_P.PHEN'${pheno}'.glm.linear' > ${GENOTYPE_TRAIN}_P_SNP.pvalue.PHEN${pheno} # create a file with SNP ID and p values
        # plink commands to build PRS score
        plink2 \
            --bfile ${GENOTYPE_TEST} \
            --score ${GENOTYPE_TRAIN}'_P.PHEN'${pheno}'.glm.linear' 3 4 12 header \
            --q-score-range range_list ${GENOTYPE_TRAIN}_P_SNP.pvalue.PHEN${pheno} \
            --out ${GENOTYPE_TEST}'_P'_PHEN${pheno}
        awk 'BEGIN { OFS="\t" } {print $3,$15}' ${GENOTYPE_TRAIN}'_MG.PHEN'${pheno}'.glm.linear' > ${GENOTYPE_TRAIN}_MG_SNP.pvalue.PHEN${pheno} # create a file with SNP ID and p values
        plink2 \
            --bfile ${GENOTYPE_TEST} \
            --score ${GENOTYPE_TRAIN}'_MG.PHEN'${pheno}'.glm.linear' 3 4 12 header \
            --q-score-range range_list ${GENOTYPE_TRAIN}_MG_SNP.pvalue.PHEN${pheno} \
            --out ${GENOTYPE_TEST}'_MG'_PHEN${pheno}
        # Append PRS results for current fold to the final files
        for thresh in "${thresholds[@]}"; do
            awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $6}' ${GENOTYPE_TEST}'_P'_PHEN${pheno}.${thresh}.sscore >> ${OUT_DIR}prs_predictions_pheno${pheno}_thresh${thresh}rep${rep}_P.txt
            awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $6}' ${GENOTYPE_TEST}'_MG'_PHEN${pheno}.${thresh}.sscore >> ${OUT_DIR}prs_predictions_pheno${pheno}_thresh${thresh}rep${rep}_MG.txt
        done
    done

    # remove intermediate files from this fold
    rm ${GENOTYPE_TRAIN}*
    rm ${GENOTYPE_TEST}*
done