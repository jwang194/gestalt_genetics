library(ggplot2)
library(data.table)
library(dplyr)
library(scales)

two_colors = c("#188085", "maroon")

polygenicity=c('0.01', '0.05', '0.1', '0.2')
shared_fraction=c('0.05', '0.2', '0.4')
genetic_covariance=c('0.01', '0.04', '0.08') # of shared components (total genetic covariance between traits is 1/2 of this value)
environmental_covariance=c('-0.1', '0.1')

num_traits=10
N=50000
M=98163
R=20

all_power_results = data.frame("polygenicity" = c(NA), "shared_fraction" = c(NA), "genetic_covariance" = c(NA), "environmental_covariance" = c(NA), "model" = c(NA), "effect_heritability" = c(NA), "power_mean" = c(NA), "power_se" = c(NA))

for (poly in polygenicity){
  for (shared in shared_fraction){
    for (cov in genetic_covariance){
      for (env in environmental_covariance){
        setting=paste0("traits_", num_traits, "_causal_", poly, "_shared_", shared, "_uniform_rg_", cov, "_random_re_", env)
        GWAS_DIR=paste0("/u/home/k/kaia/GESTALT/sim/gwas_results/", setting)
        GWAS_POWER_FILE=paste0(GWAS_DIR, "/", N, "_", M, "_", setting,"_POWER.txt" )
        df_power = fread(GWAS_POWER_FILE)
        df_power <- df_power %>%
          mutate(
            effect_heritability = ifelse(
              effect == "SPECIFIC",
              paste(effect, heritability, sep = "-"),
              "SHARED"
            )
          )
        ggplot(df_power, aes(x = model, y = power, fill = effect_heritability)) +
          geom_boxplot() +
          labs(
            x = "Model",
            y = "Power",
            fill = "Effect-Heritability"
          ) +
          scale_fill_manual(values = two_colors) + 
          theme_classic(base_size = 16)
        ggsave(filename = paste0("/u/home/k/kaia/GESTALT/sim/gwas_plots/", N, "_", M, "_", setting,"_POWER.png" ), width = 7, height = 6)
        grouped_df_power <- df_power %>%
          group_by(model, effect_heritability) %>%
          summarise(mean_power = mean(power, na.rm = TRUE), sem_power = sd(power, na.rm = TRUE) / sqrt(n()), .groups = "drop")
        grouped_df_power <- cbind(rep(poly, times = nrow(grouped_df_power)), rep(shared, times = nrow(grouped_df_power)), rep(cov, times = nrow(grouped_df_power)), rep(env, times = nrow(grouped_df_power)), grouped_df_power)
        colnames(grouped_df_power) = colnames(all_power_results)
        all_power_results <- rbind(all_power_results, grouped_df_power)
      }
    }
  }
}

all_power_results = na.omit(all_power_results)
write.table(all_power_results, file = "/u/home/k/kaia/GESTALT/sim/gwas_plots/gwas_results.txt", quote = F, row.names = F, col.names = T, sep = '\t')

# calculating per-snp effect sizes 
per_snp_variance <- function(M, polygenicity, group_fraction, heritability){
  # M = total number of variants
  # polygenicity = fraction of total variants that are causal 
  # group_fraction = fraction of the causal variants in this group (group here is shared or specific)
  # heritability = heritability of the trait 
  
  num_causal = M*polygenicity 
  num_group = num_causal*group_fraction
  persnp_variance = (heritability/num_group)*0.5 # half because effect sizes are simulated such that shared snps contribute equal heritability as specific snps (1/2 of heritability from each)
  return(persnp_variance)
}

# plots for shared variants, environmental covariance -0.1, comparing MAXH and BASE models 
M = 100000
heritability = 0.1
models_comparison = c("MAXH", "BASE")
env = -0.1
effect = "SHARED"

all_power_results_subset <- filter(all_power_results, model %in% models_comparison & effect_heritability == effect & environmental_covariance == env)
all_power_results_subset$polygenicity <- as.numeric(all_power_results_subset$polygenicity)
all_power_results_subset$shared_fraction <- as.numeric(all_power_results_subset$shared_fraction)
all_power_results_subset$genetic_covariance <- as.numeric(all_power_results_subset$genetic_covariance)
all_power_results_subset <- all_power_results_subset %>%
  mutate(persnp_variance = per_snp_variance(M, polygenicity, shared_fraction, heritability), genetic_correlation = genetic_covariance/heritability)
plot_df <- all_power_results_subset %>%
  group_by(persnp_variance, genetic_correlation) %>%
  summarize(
    power_base = mean(power_mean[model == "BASE"]), # mean because sometimes different settings (polygenicity and shared fraction) correspond to the same per-snp effect size variance, so just average those
    power_maxh = mean(power_mean[model == "MAXH"]),
    .groups = "drop"
  ) %>%
  mutate(diff = (power_maxh - power_base))

write.table(plot_df, file = paste0("/u/home/k/kaia/GESTALT/sim/gwas_plots/", N, "_", M, "_heatmap_POWER_SHARED_data.txt" ), quote = F, row.names = F, col.names = T, sep = '\t')

plot_df <- plot_df %>%
  mutate(
    # Convert to scientific notation, then order factor levels by numeric values
    persnp_variance = factor(
      format(persnp_variance, scientific = TRUE),
      levels = format(sort(unique(persnp_variance)), scientific = TRUE)
    ),
    genetic_correlation = as.factor(genetic_correlation)
  )

# Create the heatmap
ggplot(plot_df, aes(x = persnp_variance, y = genetic_correlation, fill = diff)) +
  geom_tile() +
  geom_text(aes(label = round(diff, 2)), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  #scale_fill_gradient(low = "blue", high = "red") +
  labs(
    x = "Per SNP Effect Size Variance",
    y = "Genetic Correlation of Shared Component",
    fill = "MAXH POWER - BASE POWER"
  ) + 
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0("/u/home/k/kaia/GESTALT/sim/gwas_plots/", N, "_", M, "_heatmap_POWER_SHARED.png" ), width = 10, height = 6)


# plots for shared variants, environmental covariance -0.1, comparing MAXH and BASE models 
M = 100000
heritability = 0.1
models_comparison = c("MAXH", "BASE")
env = -0.1
effect = "SPECIFIC-0.1"
num_traits = 10

all_power_results_subset <- filter(all_power_results, model %in% models_comparison & effect_heritability == effect & environmental_covariance == env)
all_power_results_subset$polygenicity <- as.numeric(all_power_results_subset$polygenicity)
all_power_results_subset$shared_fraction <- as.numeric(all_power_results_subset$shared_fraction)
all_power_results_subset$specific_fraction <- (1 - all_power_results_subset$shared_fraction)/num_traits # compute the specific fraction here
all_power_results_subset$genetic_covariance <- as.numeric(all_power_results_subset$genetic_covariance)
all_power_results_subset <- all_power_results_subset %>%
  mutate(persnp_variance = per_snp_variance(M, polygenicity, specific_fraction, heritability), genetic_correlation = genetic_covariance/heritability) # per-snp heritability for specific variants
plot_df <- all_power_results_subset %>%
  group_by(persnp_variance, genetic_correlation) %>%
  summarize(
    power_base = mean(power_mean[model == "BASE"]), # mean because sometimes different settings (polygenicity and shared fraction) correspond to the same per-snp effect size variance, so just average those
    power_maxh = mean(power_mean[model == "MAXH"]),
    .groups = "drop"
  ) %>%
  mutate(diff = (power_maxh - power_base))

write.table(plot_df, file = paste0("/u/home/k/kaia/GESTALT/sim/gwas_plots/", N, "_", M, "_heatmap_POWER_SPECIFIC_data.txt" ), quote = F, row.names = F, col.names = T, sep = '\t')


plot_df <- plot_df %>%
  mutate(
    # Convert to scientific notation, then order factor levels by numeric values
    persnp_variance = factor(
      format(persnp_variance, scientific = TRUE),
      levels = format(sort(unique(persnp_variance)), scientific = TRUE)
    ),
    genetic_correlation = as.factor(genetic_correlation)
  )

# Create the heatmap
ggplot(plot_df, aes(x = persnp_variance, y = genetic_correlation, fill = diff)) +
  geom_tile()+
  geom_text(aes(label = round(diff, 2)), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  #scale_fill_gradient(low = "blue", high = "red") +
  labs(
    x = "Per SNP Effect Size Variance",
    y = "Genetic Correlation of Shared Component",
    fill = "MAXH POWER - BASE POWER"
  ) + 
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0("/u/home/k/kaia/GESTALT/sim/gwas_plots/", N, "_", M, "_heatmap_POWER_SPECIFIC.png" ), width = 10, height = 6)



