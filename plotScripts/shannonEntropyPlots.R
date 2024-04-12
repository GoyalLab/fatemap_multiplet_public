### plotting Shannon Diversity of different datasets as a proxy for heterogeneity for Zhang, Melzer et al 2024
### Created by Madeline E Melzer on 20240210
### Last edited by Madeline E Melzer on 20240320

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(ggpubr)
library(cowplot)
set.seed(23)

resultsDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/shannonEntropy/"
plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/heterogeneity/shannonEntropy/"

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
             "TREX")

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


data = read_csv(glue("{resultsDirectory}shannonEntropy_2.csv"))

data$normDiversity = data$uniformDiversity - data$clusterDiversity

# long_data <- pivot_longer(data, 
#                           cols = c(clusterDiversity, shuffledDiversity, uniformDiversity),
#                           names_to = "DiversityType", 
#                           values_to = "DiversityValue")
# 
# long_data_filtered <- long_data %>% filter(DiversityValue > 0.51)

# Plotting Shannon Diversity Index
plot <- ggplot(data, aes(y = normDiversity, x = dataset, col = dataset)) + 
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, shape = 16) +
  #stat_compare_means(size = 4, label.x.npc = "left", paired = TRUE) +
  #facet_wrap(~resolution) +
  labs(x = "", y = "Shannon Diversity Index") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot = plot,filename=paste0(plotDirectory, "shannonDiversityIndex_resolutions.svg"), width= 18, height =12)
ggsave(plot = plot,filename=paste0(plotDirectory, "shannonDiversityIndex_resolutions.png"), width= 18, height =12)




##########################################################################################################################
# plot vs auprc
##########################################################################################################################

#data = read_csv(glue("{resultsDirectory}shannonEntropy.csv"))
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

data_summarized = data %>%
  filter(numClusters == 7) %>%
  group_by(sample, dataset) %>%
  summarise(
    clusterDiversity = mean(clusterDiversity, na.rm = TRUE),
    normDiversity = mean(normDiversity),
    uniformDiversity = mean(uniformDiversity),
    .groups = 'drop'
  )


data_plot = inner_join(data_summarized, benchmarking, by = c("sample", "dataset")) %>% mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
                                                                                             condition = str_replace(condition, "scrublet", "Scrublet"),
                                                                                             condition = str_replace(condition, "scDblFinder", "scDblFinder"),
                                                                                             condition = as.factor(condition))


ggplot(data_plot, aes(x = clusterDiversity, y = auprc, color = condition, group = condition)) +
  geom_line() + 
  theme_classic() +
  labs(x = "shannon diversity", y = "auprc", title = "auprc vs shannon diversity by detection method") +
  theme(legend.position = "bottom") 

write.csv(data_plot, glue("{resultsDirectory}ShannonEntropy_allDataForAUPRCLinePlot.csv"))








### [plotting for RANK]

scale_rank <- function(df) {
  # Rank by Phenotypic Volume
  df$rank_clusterDiversity <- rank(df$clusterDiversity, ties.method= "first")
  max_rank_clusterDiversity <- max(df$clusterDiversity)
  df$scaled_rank_clusterDiversity <- (df$clusterDiversity - 1) / (max_rank_clusterDiversity - 1) * 100
  
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

write.csv(dataRank, paste0(dataDirectory, "rankedShannonClusterDiversityForPlot.csv"))

p <- ggplot(dataRank, aes(x = rank_clusterDiversity)) +
  theme_classic() +
  labs(x = "clusterDiversity Rank", y = "AUPRC Rank", title = "")

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


ggsave(p, file = paste0(plotDirectory, 'ShannonClusterDiversityRank_AUPRCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(plotDirectory, 'ShannonClusterDiversityRank_AUPRCRank.png'), width = 5, height = 4)



# get the R squared values


lm_scDblFinder <- lm(rank_auprc_scDblFinder ~ rank_clusterDiversity, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.002014423

lm_DoubletFinder <- lm(rank_auprc_DoubletFinder ~ rank_clusterDiversity, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.004171449

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auprc_Scrublet ~ rank_clusterDiversity, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.009706617

# Linear model for hybrid
lm_hybrid <- lm(rank_auprc_hybrid ~ rank_clusterDiversity, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.01014282



##### AUROC


scale_rank_auroc <- function(df) {
  # Rank by clusterDiversity
  df$rank_clusterDiversity <- rank(df$clusterDiversity, ties.method= "first")
  max_rank_clusterDiversity <- max(df$rank_clusterDiversity)
  df$scaled_rank_clusterDiversity <- (df$rank_clusterDiversity - 1) / (max_rank_clusterDiversity - 1) * 100
  
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

write.csv(dataRank, paste0(dataDirectory, "rankedShannonClusterDiversityForPlot_AUROC.csv"))

p <- ggplot(dataRank, aes(x = rank_clusterDiversity)) +
  theme_classic() +
  labs(x = "clusterDiversity Rank", y = "AUROC Rank", title = "")

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


ggsave(p, file = paste0(plotDirectory, 'ShannonClusterDiversityRank_AUROCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(plotDirectory, 'ShannonClusterDiversityRank_AUROCRank.png'), width = 5, height = 4)



# get the R squared values


lm_scDblFinder <- lm(rank_auroc_scDblFinder ~ rank_clusterDiversity, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.005220003

lm_DoubletFinder <- lm(rank_auroc_DoubletFinder ~ rank_clusterDiversity, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.00938576

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auroc_Scrublet ~ rank_clusterDiversity, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.006830882

# Linear model for hybrid
lm_hybrid <- lm(rank_auroc_hybrid ~ rank_clusterDiversity, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.0003882647









