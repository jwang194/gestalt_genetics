# Evaluate linear combination of PRS to predict maxH phenotype
# Need to split the training fold to train GWAS effect sizes and PRS mixing weights

source ~/.bashrc

rep=$1
GENOTYPE_PREFIX="/u/home/k/kaia/GESTALT/data/sim/genotypes/N50kM100k/sim1" # Path to PLINK binary files without extension
PHENOTYPE_FILE="/u/home/k/kaia/GESTALT/data/sim/phenos/N50kM10k_10traits_rg0.08_re-0.1/" # Path to true phenotype file
NUM_FOLDS=5 # Number of folds for cross-validation
OUT_DIR="/u/home/k/kaia/GESTALT/data/sim/PRS/N50kM10k_10traits_rg0.08_re-0.1_all/" # Output directory
setting=50000_10000_traits_10_shared_0.2_specific_0.08_uniform_rg_0.08_random_re_-0.1_rep${rep} # basically prefix of the phenotype files
M=10000 # Number of snps
NUM_PHENOS=10 # Number of phenotypes
NUM_THREADS=$(nproc)

# Create sum of traits file - target for sum of PRS 
python3 sumPhenos.py ${PHENOTYPE_FILE}${setting}

# Split the fam file into N random parts (testing fold for every iteration)
mkdir ${OUT_DIR}
mkdir ${OUT_DIR}cv
shuf "${GENOTYPE_PREFIX}.fam" | split -l $(( $(wc -l < "${GENOTYPE_PREFIX}.fam") / $NUM_FOLDS )) - "${OUT_DIR}cv/fam_fold_"
num=1
for file in "${OUT_DIR}cv/fam_fold_"*; do
    mv "$file" "${OUT_DIR}cv/fam_fold${num}.IDS" # Renaming shuffled and split fam files
    ((num++))
done

# P-value thresholds for PRS scores (NOTE: in simulations all snps are independent so no LD pruning)
thresholds=(0.000005 0.001 0.05) # Add p-value thresholds here if wanted 
> range_list
for thresh in "${thresholds[@]}"; do
    echo "$thresh 0 $thresh" >> range_list # input for plink later
done

# Initialize file to store final PRS predictions (PRS predictions of testing fold will be appended too these files after every iteration)
for thresh in "${thresholds[@]}"; do
    for pheno in $(seq 1 $NUM_PHENOS); do
        echo -e "FID\tIID\tPRS" > ${OUT_DIR}prs_predictions_pheno${pheno}_thresh${thresh}rep${rep}_P.txt
    done
    echo -e "FID\tIID\tPRS" > ${OUT_DIR}prs_predictions_thresh${thresh}rep${rep}_MG.txt
    echo -e "FID\tIID\tPRS" > ${OUT_DIR}prs_predictions_thresh${thresh}rep${rep}_MG_sumPRS.txt
    echo -e "FID\tIID\tPRS" > ${OUT_DIR}prs_predictions_thresh${thresh}rep${rep}_PRS_SUM.txt
done

# Loop through each fold
for test_fold in $(seq 1 $NUM_FOLDS); do
    echo "Processing fold $test_fold..."

    GENOTYPE_TRAIN=${OUT_DIR}cv/train_fold${test_fold}
    GENOTYPE_TEST=${OUT_DIR}cv/test_fold${test_fold}

    # Create plink files for train and test
    plink2 \
     --bfile $GENOTYPE_PREFIX \
     --snps "SNP_0-SNP_"$(($M-1)) \
     --remove ${OUT_DIR}cv/fam_fold${test_fold}.IDS \
     --make-bed \
     --out ${GENOTYPE_TRAIN} \
     --threads ${NUM_THREADS}

    plink2 \
    --bfile $GENOTYPE_PREFIX \
    --snps "SNP_0-SNP_"$(($M-1)) \
    --keep ${OUT_DIR}cv/fam_fold${test_fold}.IDS \
    --make-bed \
    --out ${GENOTYPE_TEST} \
    --threads ${NUM_THREADS}

    # Filter phenotype file to keep only training individuals
    python3 ../filter_phenotype.py ${PHENOTYPE_FILE}${setting} $GENOTYPE_TRAIN 1

    # --- Call another script here to compute sum of PRS weights
    bash SumOfPRS_MaxH.sh $GENOTYPE_TRAIN $GENOTYPE_PREFIX $NUM_PHENOS "${thresholds[@]}"

    # --- GWAS on training individuals (NOTE: At this point, SNPs are already subset to the set used to simulate phenotype)
    plink2 \
    --bfile $GENOTYPE_TRAIN \
    --read-freq ${GENOTYPE_PREFIX}'.acount' \
    --pheno ${GENOTYPE_TRAIN}'_P.txt' \
    --glm allow-no-covars \
    --out ${GENOTYPE_TRAIN}'_P' \
    --threads ${NUM_THREADS}

    plink2 \
    --bfile $GENOTYPE_TRAIN \
    --read-freq ${GENOTYPE_PREFIX}'.acount' \
    --pheno ${GENOTYPE_TRAIN}'_MG.txt' \
    --glm allow-no-covars \
    --out ${GENOTYPE_TRAIN}'_MG' \
    --threads ${NUM_THREADS}

    plink2 \
    --bfile $GENOTYPE_TRAIN \
    --read-freq ${GENOTYPE_PREFIX}'.acount' \
    --pheno ${GENOTYPE_TRAIN}'_SUM.txt' \
    --glm allow-no-covars \
    --out ${GENOTYPE_TRAIN}'_SUM' \
    --no-fam-pheno \
    --threads ${NUM_THREADS}

    # ---

    for pheno in $(seq 1 $NUM_PHENOS); do
        awk 'BEGIN { OFS="\t" } {print $3,$15}' ${GENOTYPE_TRAIN}'_P.PHEN'${pheno}'.glm.linear' > ${GENOTYPE_TRAIN}_P_SNP.pvalue.PHEN${pheno} # Create a file with SNP ID and p values
        # plink commands to build PRS score for each phenotype
        plink2 \
            --bfile ${GENOTYPE_TEST} \
            --score ${GENOTYPE_TRAIN}'_P.PHEN'${pheno}'.glm.linear' 3 4 12 header \
            --q-score-range range_list ${GENOTYPE_TRAIN}_P_SNP.pvalue.PHEN${pheno} \
            --out ${GENOTYPE_TEST}'_P_PHEN'${pheno} \
            --threads ${NUM_THREADS}
        # Append PRS results for current fold to the final files
        for thresh in "${thresholds[@]}"; do
            awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $6}' ${GENOTYPE_TEST}'_P_PHEN'${pheno}.${thresh}.sscore >> ${OUT_DIR}prs_predictions_pheno${pheno}_thresh${thresh}rep${rep}_P.txt
        done
    done

    # plink command to build PRS for MaxH phenotype, only for last phenotype with is the max heritable component (first PC)
    awk 'BEGIN { OFS="\t" } {print $3,$15}' ${GENOTYPE_TRAIN}'_MG.PHEN'${NUM_PHENOS}'.glm.linear' > ${GENOTYPE_TRAIN}_MG_SNP.pvalue.PHEN${NUM_PHENOS} # Create a file with SNP ID and p values
    plink2 \
        --bfile ${GENOTYPE_TEST} \
        --score ${GENOTYPE_TRAIN}'_MG.PHEN'${NUM_PHENOS}'.glm.linear' 3 4 12 header \
        --q-score-range range_list ${GENOTYPE_TRAIN}_MG_SNP.pvalue.PHEN${NUM_PHENOS} \
        --out ${GENOTYPE_TEST}'_MG_PHEN'${NUM_PHENOS} \
        --threads ${NUM_THREADS}
    # Append PRS results for current fold to the final files
    for thresh in "${thresholds[@]}"; do
        awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $6}' ${GENOTYPE_TEST}'_MG_PHEN'${NUM_PHENOS}.${thresh}.sscore >> ${OUT_DIR}prs_predictions_thresh${thresh}rep${rep}_MG.txt
    done

    # plink command to build PRS for sum of traits (PRS OF SUM)
    awk 'BEGIN { OFS="\t" } {print $3,$15}' ${GENOTYPE_TRAIN}'_SUM.PHENSUM.glm.linear' > ${GENOTYPE_TRAIN}_SUM_SNP.pvalue.PHENSUM # Create a file with SNP ID and p values
    plink2 \
        --bfile ${GENOTYPE_TEST} \
        --score ${GENOTYPE_TRAIN}'_SUM.PHENSUM.glm.linear' 3 4 12 header \
        --q-score-range range_list ${GENOTYPE_TRAIN}_SUM_SNP.pvalue.PHENSUM \
        --out ${GENOTYPE_TEST}'_PHENSUM' \
        --threads ${NUM_THREADS}
    # Append PRS results for current fold to the final files
    for thresh in "${thresholds[@]}"; do
        awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $6}' ${GENOTYPE_TEST}'_PHENSUM'.${thresh}.sscore >> ${OUT_DIR}prs_predictions_thresh${thresh}rep${rep}_PRS_SUM.txt
    done

    # use SUM PRS betas to predict on training cohort (SUM OF PRS)
    for thresh in "${thresholds[@]}"; do
        plink2 \
            --bfile ${GENOTYPE_TEST} \
            --score ${GENOTYPE_TRAIN}_thresh_${thresh}_SUM_PRS_BETAS.txt 1 2 3 header \
            --out ${GENOTYPE_TEST}'_SUM_PRS_MG_PHEN'${NUM_PHENOS} \
            --threads ${NUM_THREADS}
        awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $6}' ${GENOTYPE_TEST}'_SUM_PRS_MG_PHEN'${NUM_PHENOS}.sscore >> ${OUT_DIR}prs_predictions_thresh${thresh}rep${rep}_MG_sumPRS.txt
    done

    # remove intermediate files from this fold
    rm -rf ${GENOTYPE_TRAIN}*
    rm -rf ${GENOTYPE_TEST}*
done

# sum the PRS scores of each phenotype to create sum of PRS, creates _SUM_sumPRS.txt
for thresh in "${thresholds[@]}"; do
        python3 sumPRS.py ${OUT_DIR} $NUM_PHENOS $thresh $rep
done


