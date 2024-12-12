library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
PRS_file_dir <- args[1] # directory to PRS predictions ex. (/u/home/k/kaia/GESTALT/data/sim/PRS/N50kM10k_5traits_rg0.1_re-0.1/)
Pheno_file_prefix <- args[2] # phenotype file dir + prefix ex. (/u/home/k/kaia/GESTALT/data/sim/phenos/N50kM10k_5traits_rg0.1_re-0.1/50000_10000_all_and_ind_overlaps_uniform_gg_random_ge_0.1_rep)
Num_pheno <- as.numeric(args[3])
Pheno_type <- args[4] # P or MG or MG_sumPRS or SUM_sumPRS or PRS_SUM
out_dir <- args[5]

if (Pheno_type == "MG_sumPRS" | Pheno_type == "MG"){
    Pheno_ref <- "MG" # define which phenotype the model tried to predict
} else if (Pheno_type == "SUM_sumPRS" | Pheno_type == "PRS_SUM") {
    Pheno_ref <- "SUM"
} else{
    Pheno_ref <- Pheno_type
}

p_threshholds <- c('0.000005', '0.001', '0.05')
R <- 25

for (thresh in p_threshholds){ # for each threshold 

    if (Pheno_type == 'P'){
        df_results <- data.frame(matrix(ncol = Num_pheno, nrow = 1)) # will append to this dataframe
        colnames(df_results) <- paste0("PRS_", 1:Num_pheno)
    }else{
        df_results <- data.frame(matrix(ncol = 1, nrow = 1)) # will append to this dataframe
        colnames(df_results) <- "PRS"
    }

    for (r in 1:R){

        if (file.exists(paste0(Pheno_file_prefix, r, "_", Pheno_ref, ".txt"))){
            df_Pheno <- as.data.frame(fread(paste0(Pheno_file_prefix, r, "_", Pheno_ref, ".txt"))) # read phenotype file

            # --- collect PRS results across phenotypes (for this threshold/iteration pair)
            if (Pheno_type == 'P'){ # if P, collect PRS predictions for each phenotype
                PRS_df <- fread(paste0(PRS_file_dir, "prs_predictions_pheno", 1, "_thresh", thresh, "rep", r, "_", Pheno_type, ".txt"))
                colnames(PRS_df)[3] <- paste0("PRS_", 1)
                for (pheno in 2:Num_pheno){
                    PRS_file <- paste0(PRS_file_dir, "prs_predictions_pheno", pheno, "_thresh", thresh, "rep", r, "_", Pheno_type, ".txt")
                    df_prs <- fread(PRS_file)
                    colnames(df_prs)[3] <- paste0("PRS_", pheno)
                    PRS_df <- merge(PRS_df, df_prs, by = c('FID', 'IID'))
                }
            }else{ # if MG or MG_sumPRS or SUM_sumPRS, just collect the PRS prediction for most heritable phenotype 
                PRS_df <- fread(paste0(PRS_file_dir, "prs_predictions_thresh", thresh, "rep", r, "_", Pheno_type, ".txt"))
                colnames(PRS_df)[3] <- "PRS"
            }
            PRS_df <- as.data.frame(PRS_df[match(df_Pheno$IID, PRS_df$IID), ])
            # ---

            # --- compute R-squared between predicted and actual
            r_squares <- c()
            if (Pheno_type == 'P'){ # if P compute R-squared between predicted and actual for each phenotype
                for (i in 3:(2 + Num_pheno)){
                    if (sum(is.na(PRS_df[,i])) == 0){
                        r_squares <- append(r_squares, summary(lm(df_Pheno[,i] ~ PRS_df[,i]))$r.squared )
                    }else{
                        r_squares <- append(r_squares, NA )
                    }
                }
            } else if (Pheno_type == 'SUM_sumPRS' | Pheno_type == "PRS_SUM") { # if SUM_sumPRS or PRS_SUM compute R-squared between predicted and actual SUM
                if (sum(is.na(PRS_df[,3])) == 0){
                    r_squares <- append(r_squares, summary(lm(df_Pheno[,3] ~ PRS_df[,3]))$r.squared )
                }else{
                    r_squares <- append(r_squares, NA )
                }
            } else{ # if MG or MG_sumPRS, compute R-squared between predicted and actual for maxh phenotype
                if (sum(is.na(PRS_df[,3])) == 0){
                    r_squares <- append(r_squares, summary(lm(df_Pheno[,(2 + Num_pheno)] ~ PRS_df[,3]))$r.squared )
                }else{
                    r_squares <- append(r_squares, NA )
                }
            }
            # --- 

            df_results <- rbind(df_results, r_squares)
        }

    }

    df_results <- df_results[-1, , drop = FALSE] # initializes dataframe with a row of NAs, so delete that here
    write.table(df_results, file = paste0(out_dir, "PRS_results_rsquared_threshold", thresh, "pheno_", Pheno_type, ".txt"), quote = F, col.names = T, row.names = F, sep = '\t')

}


