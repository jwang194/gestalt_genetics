library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)


setting <- "C:/Users/kaiak/OneDrive/Documents/GESTALT/N50kM10k_10traits_rg0.05_reneg0.1/"
p_thresholds <- c('0.000005', '0.001', '0.05')
num_pheno <- 10
#heritabilities <- c(0.4, 0.3, 0.15, 0.5, 0.2, 'MaxH')
heritabilities <- c(rep(0.1, times = num_pheno), 'MaxH')

# --- plot all p value settings separately
for (p_current in p_thresholds){
  df_P <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_P.txt"))
  df_MG <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_MG.txt"), select = c(num_pheno))
  colnames(df_MG) <- "PRS_MaxH"
  df_merged <- data.frame(df_P, df_MG)
  long_df <- stack(df_merged[, 1:(num_pheno+1)])
  colnames(long_df) <- c("Rsquared", "Trait")
  long_df$Heritability <- rep(heritabilities, each = nrow(df_merged))
  
  ggplot(long_df, aes(x = Heritability, y = Rsquared)) + 
    geom_boxplot(color = "darkblue")+
    labs(x = "Heritability of Trait", y = expression(paste(R^2))) + 
    theme_bw() +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(setting, "threshold_" ,p_current , "_PRS.png"), width = 6, height = 6)
}


# --- plot best p value cutoff for each trait
df_P <- fread(paste0(setting, "PRS_results_rsquared_threshold", '0.000005', "pheno_P.txt"))
df_MG <- fread(paste0(setting, "PRS_results_rsquared_threshold", '0.000005', "pheno_MG.txt"), select = c(num_pheno))
colnames(df_MG) <- "PRS_MaxH"
df_best <- data.frame(df_P, df_MG)
for (p_current in p_thresholds[2:length(p_thresholds)]){
  df_P <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_P.txt"))
  df_MG <- fread(paste0(setting, "PRS_results_rsquared_threshold", p_current, "pheno_MG.txt"), select = c(num_pheno))
  colnames(df_MG) <- "PRS_MaxH"
  df_merged <- data.frame(df_P, df_MG)
  for (i in 1:ncol(df_merged)){
    if (mean(df_merged[,i]) > mean(df_best[,i])){
      df_best[,i] <- df_merged[,i]
    }
  }
}
long_df <- stack(df_best[, 1:(num_pheno+1)])
colnames(long_df) <- c("Rsquared", "Trait")
long_df$Heritability <- rep(heritabilities, each = nrow(df_best))

ggplot(long_df, aes(x = Heritability, y = Rsquared)) + 
  geom_boxplot(color = "darkblue")+
  labs(x = "Heritability of Trait", y = expression(paste(R^2))) + 
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = paste0(setting, "threshold_best_PRS.png"), width = 6, height = 6)
