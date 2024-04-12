### plotting differential expresssion results for heterogeneity of samples for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20240215
### Last edited by Madeline E Melzer on 20240320

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(PRROC)
library(conflicted)


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


DEDataDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/DE/data"
dePlotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/heterogeneity/DEHeterogeneity/"

data= read.csv(glue("{DEDataDir}/DEResults_wilcox.csv"))


# this is just to add the number of clusters
data <- data %>%
  group_by(Dataset, Sample, Resolution) %>%
  mutate(numClusters = n_distinct(Cluster)) %>%
  ungroup()


diff_data <- data %>%
  filter(numClusters == 7) %>%
  group_by(Dataset, Sample) %>%  # Group by dataset and sample or any other relevant variable
  summarize(pDiff = mean(pUp + pDown, na.rm = TRUE),
            fDiff = mean(fU - fL, na.rm = TRUE),
            .groups = 'drop')  # Compute the mean difference pUp - pDown

plot = ggplot(diff_data, aes(x = Dataset, y = pDiff, fill = Dataset)) +
  geom_boxplot(aes(fill = Dataset), color = "black", position = position_dodge(width = 0.9), width = 0.85) +
  #geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(x = "Dataset",
       y = "Average pUp + pDown") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
plot

ggsave(plot, file = paste0(dePlotDirectory, 'pUpPluspDownPerDataset_0.6res.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(dePlotDirectory, 'pUpPluspDownPerDataset_0.6res.png'), width = 6, height = 4)



plot = ggplot(diff_data, aes(x = Dataset, y = fDiff, fill = Dataset)) +
  geom_boxplot(aes(fill = Dataset), color = "black", position = position_dodge(width = 0.9), width = 0.85) +
  #geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(x = "Dataset",
       y = "Average fU - fL Difference") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
plot

ggsave(plot, file = paste0(dePlotDirectory, 'upperLowerBounds_0.6res.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(dePlotDirectory, 'upperLowerBounds_0.6res.png'), width = 6, height = 4)



####### plotting AUPRC vs DE metrics
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


diff_data = diff_data %>% dplyr::rename(dataset = Dataset, sample = Sample)

data_plot = inner_join(diff_data, benchmarking, by = c("sample", "dataset")) %>% mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
                                                                                        condition = str_replace(condition, "scrublet", "Scrublet"),
                                                                                        condition = str_replace(condition, "scDblFinder", "scDblFinder"),
                                                                                        condition = as.factor(condition))

ggplot(data_plot, aes(x = fDiff, y = auprc, color = condition, group = condition)) +
  geom_line() + 
  theme_classic() +
  labs(x = "Average fU - fL", y = "auprc", title = "auprc vs fDiff by detection method") +
  theme(legend.position = "bottom") 

write.csv(data_plot, glue("{DEDataDir}/DEHeterogeneity_allDataForAUPRCLinePlot.csv"))















### [plotting for RANK]

scale_rank <- function(df) {
  # Rank by DE
  df$rank_pDiff <- rank(df$pDiff, ties.method= "first")
  max_rank_pDiff <- max(df$rank_pDiff)
  df$scaled_rank_pDiff <- (df$rank_pDiff - 1) / (max_rank_pDiff - 1) * 100
  
  df$rank_fDiff <- rank(df$fDiff, ties.method= "first")
  max_rank_fDiff <- max(df$rank_fDiff)
  df$scaled_rank_fDiff <- (df$rank_fDiff - 1) / (max_rank_fDiff - 1) * 100
  
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

write.csv(dataRank, paste0(DEDataDir, "rankedDEHeterogeneityForPlot.csv"))


## pDiff
p <- ggplot(dataRank, aes(x = rank_pDiff)) +
  theme_classic() +
  labs(x = "pDiff Rank", y = "AUPRC Rank", title = "")

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


ggsave(p, file = paste0(dePlotDirectory, 'pDiffRank_AUPRCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(dePlotDirectory, 'pDiffRank_AUPRCRank.png'), width = 5, height = 4)

# get the R squared values

lm_scDblFinder <- lm(rank_auprc_scDblFinder ~ rank_pDiff, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.283045

lm_DoubletFinder <- lm(rank_auprc_DoubletFinder ~ rank_pDiff, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.002703772

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auprc_Scrublet ~ rank_pDiff, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.09836309

# Linear model for hybrid
lm_hybrid <- lm(rank_auprc_hybrid ~ rank_pDiff, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.03342072



## fDiff
p <- ggplot(dataRank, aes(x = rank_fDiff)) +
  theme_classic() +
  labs(x = "fDiff Rank", y = "AUPRC Rank", title = "")

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


ggsave(p, file = paste0(dePlotDirectory, 'fDiffRank_AUPRCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(dePlotDirectory, 'fDiffRank_AUPRCRank.png'), width = 5, height = 4)



# get the R squared values


lm_scDblFinder <- lm(rank_auprc_scDblFinder ~ rank_fDiff, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.0003882647

lm_DoubletFinder <- lm(rank_auprc_DoubletFinder ~ rank_fDiff, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.1261863

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auprc_Scrublet ~ rank_fDiff, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.1053495

# Linear model for hybrid
lm_hybrid <- lm(rank_auprc_hybrid ~ rank_fDiff, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.01654469



##### AUROC


scale_rank_auroc <- function(df) {
  # Rank by DE
  df$rank_pDiff <- rank(df$pDiff, ties.method= "first")
  max_rank_pDiff <- max(df$rank_pDiff)
  df$scaled_rank_pDiff <- (df$rank_pDiff - 1) / (max_rank_pDiff - 1) * 100
  
  df$rank_fDiff <- rank(df$fDiff, ties.method= "first")
  max_rank_fDiff <- max(df$rank_fDiff)
  df$scaled_rank_fDiff <- (df$rank_fDiff - 1) / (max_rank_fDiff - 1) * 100
  
  # Rank by AUROC
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

write.csv(dataRank, paste0(DEDataDir, "rankedDEHeterogeneityForPlot_AUROC.csv"))


## pDiff

p <- ggplot(dataRank, aes(x = rank_pDiff)) +
  theme_classic() +
  labs(x = "pDiff Rank", y = "AUROC Rank", title = "")

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


ggsave(p, file = paste0(dePlotDirectory, 'pDiffRank_AUROCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(dePlotDirectory, 'pDiffRank_AUROCRank.png'), width = 5, height = 4)


# get the R squared values

lm_scDblFinder <- lm(rank_auroc_scDblFinder ~ rank_pDiff, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.1879201

lm_DoubletFinder <- lm(rank_auroc_DoubletFinder ~ rank_pDiff, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.01948064

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auroc_Scrublet ~ rank_pDiff, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.09029581

# Linear model for hybrid
lm_hybrid <- lm(rank_auroc_hybrid ~ rank_pDiff, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.07670505





## fDiff

p <- ggplot(dataRank, aes(x = rank_fDiff)) +
  theme_classic() +
  labs(x = "fDiff Rank", y = "AUROC Rank", title = "")

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

ggsave(p, file = paste0(dePlotDirectory, 'fDiffRank_AUROCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(dePlotDirectory, 'fDiffRank_AUROCRank.png'), width = 5, height = 4)


# get the R squared values

lm_scDblFinder <- lm(rank_auroc_scDblFinder ~ rank_fDiff, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.05036057

lm_DoubletFinder <- lm(rank_auroc_DoubletFinder ~ rank_fDiff, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.1757844

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auroc_Scrublet ~ rank_fDiff, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.06845083

# Linear model for hybrid
lm_hybrid <- lm(rank_auroc_hybrid ~ rank_fDiff, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.0003463225





