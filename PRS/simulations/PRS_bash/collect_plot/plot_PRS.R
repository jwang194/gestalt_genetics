library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)


setting <- "C:/Users/kaiak/OneDrive/Documents/GESTALT/N50kM10k_10traits_rg0.05_reneg0.1/"
#setting <- "C:/Users/kaiak/OneDrive/Documents/GESTALT/N50kM10k_5traits_rg0.1_reneg0.1/"
p_thresholds <- c('0.000005', '0.001', '0.05')
num_pheno <- 10
#num_pheno <- 5
#heritabilities <- c(0.4, 0.3, 0.15, 0.5, 0.2, 'MaxH', 'SUM_MaxH', 'SUM')
heritabilities <- c(rep(0.1, times = num_pheno), 'MaxH', 'SUM_MaxH', 'SUM')

# --- plot all p value settings separately
for (p_current in p_thresholds){
  df_P <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_P.txt"))
  df_MG <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_MG.txt"))
  df_MG_sum <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_MG_sumPRS.txt"))
  df_sum_sum <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_SUM_sumPRS.txt"))
  colnames(df_MG) <- "PRS_MaxH"
  colnames(df_MG_sum) <- "PRS_SUM_MaxH"
  colnames(df_sum_sum) <- "PRS_SUM_SUM"
  df_merged <- data.frame(df_P, df_MG, df_MG_sum, df_sum_sum)
  long_df <- stack(df_merged[, 1:(num_pheno+3)])
  colnames(long_df) <- c("Rsquared", "Trait")
  long_df$Heritability <- rep(heritabilities, each = nrow(df_merged))
  
  ggplot(long_df, aes(x = Heritability, y = Rsquared)) + 
    geom_boxplot(color = "darkblue")+
    labs(x = "Heritability of Trait", y = expression(paste(R^2))) + 
    theme_bw() +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(setting, "threshold_" ,p_current , "_PRS_with_sumPRS_MaxH_and_SumTraits.png"), width = 6, height = 6)
}


# --- plot best p value cutoff for each trait
df_P <- fread(paste0(setting, "PRS_results_rsquared_threshold", '0.000005', "pheno_P.txt"))
df_MG <- fread(paste0(setting, "PRS_results_rsquared_threshold", '0.000005', "pheno_MG.txt"))
df_MG_sum <- fread(paste0(setting, "PRS_results_rsquared_threshold", '0.000005', "pheno_MG_sumPRS.txt"))
df_sum_sum <- fread(paste0(setting, "PRS_results_rsquared_threshold", '0.000005', "pheno_SUM_sumPRS.txt"))
colnames(df_MG) <- "PRS_MaxH"
colnames(df_MG_sum) <- "PRS_SUM"
colnames(df_sum_sum) <- "PRS_SUM_SUM"
df_best <- data.frame(df_P, df_MG, df_MG_sum, df_sum_sum)
for (p_current in p_thresholds[2:length(p_thresholds)]){
  df_P <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_P.txt"))
  df_MG <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_MG.txt"))
  df_MG_sum <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_MG_sumPRS.txt"))
  df_sum_sum <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_SUM_sumPRS.txt"))
  colnames(df_MG) <- "PRS_MaxH"
  colnames(df_MG_sum) <- "PRS_SUM"
  colnames(df_sum_sum) <- "PRS_SUM_SUM"
  df_merged <- data.frame(df_P, df_MG, df_MG_sum, df_sum_sum)
  for (i in 1:ncol(df_merged)){
    if (mean(df_merged[,i], na.rm = T) > mean(df_best[,i], na.rm = T)){
      df_best[,i] <- df_merged[,i]
    }
  }
}
long_df <- stack(df_best[, 1:(num_pheno+3)])
colnames(long_df) <- c("Rsquared", "Trait")
long_df$Heritability <- rep(heritabilities, each = nrow(df_best))

ggplot(long_df, aes(x = Heritability, y = Rsquared)) + 
  geom_boxplot(color = "darkblue")+
  labs(x = "Heritability of Trait", y = expression(paste(R^2))) + 
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = paste0(setting, "threshold_best_PRS_with_sumPRS_MaxH_and_SumTraits.png"), width = 6, height = 6)
