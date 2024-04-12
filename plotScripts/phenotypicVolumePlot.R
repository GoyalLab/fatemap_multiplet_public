### plotting phenotypic volume of different datasets as a proxy for heterogeneity for Zhang, Melzer et al 2024
### Created by Madeline E Melzer on 20240210
### Last edited by Madeline E Melzer on 20240210

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(ggpubr)
library(cowplot)


set.seed(23)

dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/heterogeneity/phenotypicVolume/data_variableGenes/"
plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/heterogeneity/phenotypicVolume/"

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
                   "cellTag" = "CellTag-Multi",
                   "watermelon" = "Watermelon")


##########################################################################################################################
# Combining all files
##########################################################################################################################

files <- list.files(dataDirectory, pattern = "\\.csv$", full.names = TRUE)
data <- files %>%
  map_dfr(~{
    df <- read_csv(.x)
    
    filename <- basename(.x)
    parts <- strsplit(filename, "__|\\.csv")[[1]]
    dataset <- parts[1]
    sample <- parts[2]
    
    df <- df %>%
      mutate(dataset = dataset, sample = sample)
    
    return(df)
  })

write.csv(data, paste0(dataDirectory, "allSamples_PV_variableGenes.csv"))
data = read.csv(paste0(dataDirectory, "allSamples_PV_variableGenes.csv"))

data$dataset <- factor(data$dataset, levels = barcoded)
data <- data %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

plot = ggplot(data, aes(x=dataset, y=allPV)) +
  #geom_violin() +
  geom_boxplot(width=0.1) +
  stat_compare_means() + #this would give mann whitney test result
  theme_classic() +
  ylab("log of phenotypic volume(gene-gene)") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot






##########################################################################################################################
# plot vs auprc
##########################################################################################################################

data = read_csv(paste0(dataDirectory, "allSamples_PV_variableGenes.csv"))
benchmarking = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/barcodedAllForPlot.csv")

benchmarking = benchmarking %>% 
  filter(dbl_act == 0.08) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename)) %>%
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
  group_by(sample, dataset) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename)) %>%
  summarise(
    PV = mean(allPV, na.rm = TRUE),
    .groups = 'drop'      
  )

data_plot = inner_join(data_summarized, benchmarking, by = c("sample", "dataset")) %>% mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
                                                                                              condition = str_replace(condition, "scrublet", "Scrublet"),
                                                                                              condition = str_replace(condition, "scDblFinder", "scDblFinder"),
                                                                                              condition = as.factor(condition))

data_plot = data_plot %>% 
  group_by(sample, condition) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),
    PV = mean(PV, na.rm = TRUE),
    .groups = 'drop'      
  )

ggplot(data_plot, aes(x = PV, y = auprc, color = condition, group = condition)) +
  geom_line() + 
  theme_classic() +
  labs(x = "PV", y = "auprc", title = "auprc vs PV by detection method") +
  theme(legend.position = "bottom") 










### [plotting for RANK]

scale_rank <- function(df) {
  # Rank by Phenotypic Volume
  df$rank_PV <- rank(df$PV, ties.method= "first")
  max_rank_PV <- max(df$rank_PV)
  df$scaled_rank_PV <- (df$rank_PV - 1) / (max_rank_PV - 1) * 100
  
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

write.csv(dataRank, paste0(dataDirectory, "rankedPhenotypicVolumeForPlot.csv"))

p <- ggplot(dataRank, aes(x = rank_E_distance)) +
  theme_classic() +
  labs(x = "PV Rank", y = "AUPRC Rank", title = "")

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


ggsave(p, file = paste0(plotDirectory, 'phenotypicVolumeRank_AUPRCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(plotDirectory, 'phenotypicVolumeRank_AUPRCRank.png'), width = 5, height = 4)



# get the R squared values


lm_scDblFinder <- lm(rank_auprc_scDblFinder ~ rank_PV, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.0397075

lm_DoubletFinder <- lm(rank_auprc_DoubletFinder ~ rank_PV, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 0.01730478

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auprc_Scrublet ~ rank_PV, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.07880565

# Linear model for hybrid
lm_hybrid <- lm(rank_auprc_hybrid ~ rank_PV, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.0002027755



##### AUROC


scale_rank_auroc <- function(df) {
  # Rank by Phenotypic Volume
  df$rank_PV <- rank(df$PV, ties.method= "first")
  max_rank_PV <- max(df$rank_PV)
  df$scaled_rank_PV <- (df$rank_PV - 1) / (max_rank_PV - 1) * 100
  
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

write.csv(dataRank, paste0(dataDirectory, "rankedPhenotypicVolumeForPlot_AUROC.csv"))

p <- ggplot(dataRank, aes(x = rank_PV)) +
  theme_classic() +
  labs(x = "PV Rank", y = "AUROC Rank", title = "")

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


ggsave(p, file = paste0(plotDirectory, 'phenotypicVolumeRank_AUROCRank.svg'), width = 5, height = 4)
ggsave(p, file = paste0(plotDirectory, 'phenotypicVolumeRank_AUROCRank.png'), width = 5, height = 4)



# get the R squared values


lm_scDblFinder <- lm(rank_auroc_scDblFinder ~ rank_PV, data = dataRank)
r2_scDblFinder <- summary(lm_scDblFinder)$r.squared
r2_scDblFinder # 0.01884334

lm_DoubletFinder <- lm(rank_auroc_DoubletFinder ~ rank_PV, data = dataRank)
r2_DoubletFinder <- summary(lm_DoubletFinder)$r.squared
r2_DoubletFinder # 7.728535e-05

# Linear model for Scrublet
lm_Scrublet <- lm(rank_auroc_Scrublet ~ rank_PV, data = dataRank)
r2_Scrublet <- summary(lm_Scrublet)$r.squared
r2_Scrublet # 0.08972548

# Linear model for hybrid
lm_hybrid <- lm(rank_auroc_hybrid ~ rank_PV, data = dataRank)
r2_hybrid <- summary(lm_hybrid)$r.squared
r2_hybrid # 0.001117225
