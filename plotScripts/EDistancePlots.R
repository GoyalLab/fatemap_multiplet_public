### plotting E-distances of different datasets as a proxy for heterogeneity for Zhang, Melzer et al 2024
### Created by Madeline E Melzer on 20240210
### Last edited by Madeline E Melzer on 20240319

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(ggpubr)
library(cowplot)
set.seed(23)

resultsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/EDistance/results/"
plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/heterogeneity/EDistance/"

barcoded = c("FM01", 
             "FM02", 
             "FM03", 
             "FM04", 
             "FM05", 
             "FM06", 
             "FM08", 
             "Biorxiv",
             "non_cancer", 
             "ClonMapper", 
             "SPLINTR", 
             "LARRY",
             "cellTag",
             "watermelon",
             "TREX",
             "smartseq3_reads",
             "smartseq3_umis")

barcodedRename = c("FM01" = "Goyal et al. 1", 
                   "FM02" = "Goyal et al. 2", 
                   "FM03" = "Goyal et al. 3", 
                   "FM04" = "Goyal et al. 4", 
                   "FM05" = "Goyal et al. 5", 
                   "FM06" = "Goyal et al. 6", 
                   "FM08" = "Goyal et al. 8", 
                   "non_cancer" = "Jiang et al.", 
                   #"ClonMapper", 
                   "Biorxiv" = "Jain et al.", 
                   #"SPLINTR", 
                   #"TREX",
                   #"LARRY",
                   "cellTag" = "CellTag-multi",
                   "watermelon" = "Watermelon")

data = read_csv(glue("{resultsDirectory}EDistances_allDataAndControls.csv"))
data$numClusters <- as.factor(data$numClusters)

# plot all of the cluster numbers per dataset on the same axis in a stacked bar plot
ggplot(data, aes(x = dataset, y = E_distance, fill = numClusters, group = interaction(numClusters, dataset))) +
  geom_col(position = "stack") +
  theme_minimal() +
  labs(x = "Dataset", y = "E-Distance", fill = "Number of Clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_data <- data %>%
  filter(identity == "sample") %>%
  filter(!is.na(clusterA), !is.na(clusterB), clusterA != clusterB) %>%
  group_by(sample, dataset, numClusters) %>%
  summarize(avg_E_distance = mean(E_distance, na.rm = TRUE), .groups = 'drop')

common_clusters <- sample_data %>%
  dplyr::count(numClusters) %>%
  filter(n == max(n)) %>%
  pull(numClusters)

final_data <- sample_data %>%
  #filter(numClusters %in% common_clusters) %>%
  filter(numClusters == 4) %>%
  group_by(dataset) %>%
  summarize(avg_E_distance = mean(avg_E_distance, na.rm = TRUE), .groups = 'drop')

# Plot
ggplot(final_data, aes(x = dataset, y = avg_E_distance, fill = dataset)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Average EDistance by Dataset",
       x = "Dataset",
       y = "Average EDistance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### plot with controls

agg_data <- data %>%
  filter(!is.na(clusterA), !is.na(clusterB), clusterA != clusterB) %>%
  group_by(identity, dataset, numClusters, sample) %>%
  summarize(avg_E_distance = mean(E_distance, na.rm = TRUE), .groups = 'drop')


agg_data_filtered <- agg_data %>%
  filter(numClusters == 7)

ggplot(agg_data_filtered, aes(x = dataset, y = avg_E_distance, fill = identity)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Average EDistance by Dataset and Identity",
       x = "Dataset",
       y = "Average EDistance") +
  scale_fill_manual(values = c("control" = "blue", "sample" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
































####### plotting AUPRC vs E-distance
benchmarking = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/barcodedAllForPlot.csv")

benchmarking = benchmarking %>% 
  filter(dbl_act == 0.08) %>%
  filter(
    !(
      (dataset == "SPLINTR" & sample %in% c("inVitro_KRAS", "retransplant")) |
        (dataset == "TREX" & sample == "brain1") |
        (dataset == "smartseq3" & sample %in% c("sample1"))
    )
  ) %>%
  group_by(dataset, condition, sample) %>%
  mutate(sample = case_when(
    sample %in% c("1_DMSO_A", "2_DMSO_B") ~ "mutated",
    TRUE ~ sample
  )) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp and dbl_act
    .groups = 'drop'      
  )

data_filtered = data %>%
  filter(!is.na(clusterA), !is.na(clusterB), clusterA != clusterB) %>%
  filter(identity == "sample") %>%
  group_by(sample, numClusters, dataset) %>%
  summarise(
    E_distance = mean(E_distance, na.rm = TRUE),
    .groups = 'drop'      
  )

data_plot = inner_join(data_filtered, benchmarking, by = c("sample", "dataset")) %>% mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
                                                                                            condition = str_replace(condition, "scrublet", "Scrublet"),
                                                                                            condition = str_replace(condition, "scDblFinder", "scDblFinder"),
                                                                                            condition = as.factor(condition))

data_plot = data_plot %>% 
  filter(numClusters == 7) %>% # 7 and 4 are most common clusters.
  group_by(sample, condition, dataset) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),
    E_distance = mean(E_distance, na.rm = TRUE),
    .groups = 'drop'      
  )

write.csv(data_plot, glue("{resultsDirectory}EDistances_allDataForAUPRCLinePlot.csv"))



ggplot(data_plot, aes(x = E_distance, y = auprc, color = condition, group = condition)) +
  geom_line() + 
  theme_classic() +
  labs(x = "E-distance", y = "auprc", title = "auprc vs E-distance by detection method") +
  theme(legend.position = "bottom") 





### [plotting for RANK]

#data_plot = read.csv(glue("{resultsDirectory}EDistances_allDataForAUPRCLinePlot.csv"))

scale_rank <- function(df) {
  # Rank by E_distance
  df$rank_E_distance <- rank(df$E_distance, ties.method= "first")
  max_rank_E_distance <- max(df$rank_E_distance)
  df$scaled_rank_E_distance <- (df$rank_E_distance - 1) / (max_rank_E_distance - 1) * 100
  
  # Rank by AUPRC
  for (cond in unique(df$condition)) {
    # Dynamically generate new column names for ranks and scaled ranks
    rank_col <- sym(paste0("rank_auprc_", cond))
    max_rank_col <- sym(paste0("max_rank_auprc_", cond))
    scaled_rank_col <- sym(paste0("scaled_rank_auprc_", cond))
    
    df <- df %>%
      group_by(condition) %>%
      mutate(
        !!rank_col := if_else(condition == cond, rank(auprc, ties.method = "first"), NA_real_)
      ) %>%
      ungroup() %>%
      mutate(
        !!max_rank_col := max(!!rank_col, na.rm = TRUE),
        !!scaled_rank_col := if_else(!is.na(!!rank_col), ((!!rank_col - 1) / (!!max_rank_col - 1) * 100), NA_real_)
      )
  }
  # remove the intermediate max_rank columns if they are not needed
  df <- select(df, -starts_with("max_rank_auprc_"))
  
  return(df)
}

dataRank <- scale_rank(data_plot)

write.csv(dataRank, paste0(resultsDirectory, "rankedEDistanceForPlot.csv"))

p <- ggplot(dataRank, aes(x = rank_E_distance)) +
  theme_classic() +
  labs(x = "E_distance Rank", y = "AUPRC Rank", title = "")

# Manually add points for each condition
p <- p + 
  geom_point(aes(y = rank_auprc_scDblFinder, color = "scDblFinder")) +
  geom_point(aes(y = rank_auprc_DoubletFinder, color = "DoubletFinder")) +
  geom_point(aes(y = rank_auprc_Scrublet, color = "Scrublet")) +
  geom_point(aes(y = rank_auprc_hybrid, color = "hybrid"))

p <- p +
  geom_smooth(aes(y = rank_auprc_scDblFinder, color = "scDblFinder"), method = "lm", formula = y ~ x, se = FALSE) +
  geom_smooth(aes(y = rank_auprc_DoubletFinder, color = "DoubletFinder"), method = "lm", formula = y ~ x, se = FALSE) +
  geom_smooth(aes(y = rank_auprc_Scrublet, color = "Scrublet"), method = "lm", formula = y ~ x, se = FALSE) +
  geom_smooth(aes(y = rank_auprc_hybrid, color = "hybrid"), method = "lm", formula = y ~ x, se = FALSE)

# Use your custom palette
p <- p + scale_color_manual(values = palette)

# Final adjustments
p <- p + theme(legend.title = element_blank(), legend.text = element_text(size = 10))

# Display the plot
print(p)


ggsave(p, file = paste0(plotDirectory, 'EDistanceRank_AUPRCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(plotDirectory, 'EDistanceRank_AUPRCRank.png'), width = 5, height = 4)



# get the R squared values


lm_scDblFinder <- lm(rank_auprc_scDblFinder ~ rank_E_distance, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.1501722

lm_DoubletFinder <- lm(rank_auprc_DoubletFinder ~ rank_E_distance, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.006740707

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auprc_Scrublet ~ rank_E_distance, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.1301049

# Linear model for hybrid
lm_hybrid <- lm(rank_auprc_hybrid ~ rank_E_distance, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.09974119



##### AUROC


scale_rank_auroc <- function(df) {
  # Rank by E_distance
  df$rank_E_distance <- rank(df$E_distance, ties.method= "first")
  max_rank_E_distance <- max(df$rank_E_distance)
  df$scaled_rank_E_distance <- (df$rank_E_distance - 1) / (max_rank_E_distance - 1) * 100
  
  # Rank by AUPRC
  for (cond in unique(df$condition)) {
    # Dynamically generate new column names for ranks and scaled ranks
    rank_col <- sym(paste0("rank_auroc_", cond))
    max_rank_col <- sym(paste0("max_rank_auroc_", cond))
    scaled_rank_col <- sym(paste0("scaled_rank_auroc_", cond))
    
    df <- df %>%
      group_by(condition) %>%
      mutate(
        !!rank_col := if_else(condition == cond, rank(auroc, ties.method = "first"), NA_real_)
      ) %>%
      ungroup() %>%
      mutate(
        !!max_rank_col := max(!!rank_col, na.rm = TRUE),
        !!scaled_rank_col := if_else(!is.na(!!rank_col), ((!!rank_col - 1) / (!!max_rank_col - 1) * 100), NA_real_)
      )
  }
  # remove the intermediate max_rank columns if they are not needed
  df <- select(df, -starts_with("max_rank_auroc_"))
  
  return(df)
}

dataRank <- scale_rank_auroc(data_plot)

write.csv(dataRank, paste0(dataDirectory, "rankedEDistanceForPlot_AUROC.csv"))

p <- ggplot(dataRank, aes(x = rank_E_distance)) +
  theme_classic() +
  labs(x = "E_distance Rank", y = "AUROC Rank", title = "")

# Manually add points for each condition
p <- p + 
  geom_point(aes(y = rank_auroc_scDblFinder, color = "scDblFinder")) +
  geom_point(aes(y = rank_auroc_DoubletFinder, color = "DoubletFinder")) +
  geom_point(aes(y = rank_auroc_Scrublet, color = "Scrublet")) +
  geom_point(aes(y = rank_auroc_hybrid, color = "hybrid"))

p <- p +
  geom_smooth(aes(y = rank_auroc_scDblFinder, color = "scDblFinder"), method = "lm", formula = y ~ x, se = FALSE) +
  geom_smooth(aes(y = rank_auroc_DoubletFinder, color = "DoubletFinder"), method = "lm", formula = y ~ x, se = FALSE) +
  geom_smooth(aes(y = rank_auroc_Scrublet, color = "Scrublet"), method = "lm", formula = y ~ x, se = FALSE) +
  geom_smooth(aes(y = rank_auroc_hybrid, color = "hybrid"), method = "lm", formula = y ~ x, se = FALSE)

# Use your custom palette
p <- p + scale_color_manual(values = palette)

# Final adjustments
p <- p + theme(legend.title = element_blank(), legend.text = element_text(size = 10))

# Display the plot
print(p)


ggsave(p, file = paste0(plotDirectory, 'EDistanceRank_AUROCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(plotDirectory, 'EDistanceRank_AUROCRank.png'), width = 5, height = 4)



# get the R squared values


lm_scDblFinder <- lm(rank_auroc_scDblFinder ~ rank_E_distance, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.06845083

lm_DoubletFinder <- lm(rank_auroc_DoubletFinder ~ rank_E_distance, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.0008735956

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auroc_Scrublet ~ rank_E_distance, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.09836309

# Linear model for hybrid
lm_hybrid <- lm(rank_auroc_hybrid ~ rank_E_distance, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.09563565


