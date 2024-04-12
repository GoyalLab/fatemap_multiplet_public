### plotting alternate technologies (smartseq3) together for Zhang, Melzer et al 2024
### Created by Madeline E Melzer on 20240210
### Last edited by Madeline E Melzer on 20240314

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(PRROC)

plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/seqTechComparison/"
dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/"


palette = c("classifier" = "#f59fff", 
            "scDblFinder" = "#b26c98", 
            "DoubletFinder" = "#feb325", 
            "Scrublet" = "#2474a9", 
            "hybrid" = "#03bfc4", 
            "scrambled" = "#007024",
            "otherMethods" = "#808285")

################################
#                              #
#        Data Formatting       #
#                              #
################################

subfolders <- list.dirs(dataDirectory, recursive = FALSE)

benchmarkingData <- map_dfr(subfolders, ~{
  detection_file_path <- file.path(.x, "stats/all_detection_rates.tsv")
  
  if (file.exists(detection_file_path)) {
    read_tsv(detection_file_path) %>%
      separate(ID, into = c("dataset", "sample"), sep = "_", extra = "merge") %>%
      mutate(#sample = ifelse(is.na(sample), NA_character_, sample), 
        condition = basename(.x),
        auroc = roc,
        auprc = pr) %>%
      dplyr::select(dataset, sample, condition, auroc, auprc, dbl_exp, dbl_act)
  } else {
    tibble()  # Return an empty tibble if the file doesn't exist
  }
})

benchmarkingData
write.csv(benchmarkingData, paste0(dataDirectory, "allBenchmarking.csv"))
#benchmarkingData = read_csv(paste0(dataDirectory, "allBenchmarking.csv"))




benchmarkingData = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/barcodedAllForPlot.csv") #20240314

### all barcoded datasets
benchmarkingData_formatted <- benchmarkingData %>%
  filter(dataset %in% c("smartseq3_reads", "smartseq3_umis", "TREX"),
         sample != "sample1") %>%
  # mutate(dataset = case_when( #commented out this part on 20240314 when barcodedAllForPlot.csv was being used because this specification was already made. 
  #   str_detect(sample, "reads") ~ "smartseq3_reads",
  #   str_detect(sample, "umis") ~ "smartseq3_umis",
  #   TRUE ~ dataset  # Keep original dataset value if no match
  # )) %>%
  #filter(dataset != "smartseq3_umis") %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition)) %>%
  filter(dbl_act == 0.08) #%>%
  # group_by(dataset, condition, sample) %>%
  # summarise(
  #   auprc = mean(auprc, na.rm = TRUE),
  #   auroc = mean(auroc, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp and dbl_act
  #   .groups = 'drop'
  # )






  

### AUPRC
benchmarkingData_summary = benchmarkingData_formatted %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
benchmarkingData_plot <- left_join(benchmarkingData_formatted, benchmarkingData_summary, by = c("dataset", "condition"))

benchmarkingData_plot = benchmarkingData_plot %>% filter(dataset != "smartseq3_umis")

plot = ggplot(benchmarkingData_plot, aes(x = dataset, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  #geom_jitter(data = benchmarkingData_plot, aes(x = dataset, y = auprc, group = condition, color = condition), position = position_dodge(width = 0.75)) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_smartseq3_AUPRC_reads_0.08.svg'), width = 2.5, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_smartseq3_AUPRC_reads_0.08.png'), width = 2.5, height = 3)


### AUROC
benchmarkingData_summary = benchmarkingData_formatted %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc))
benchmarkingData_plot <- left_join(benchmarkingData_formatted, benchmarkingData_summary, by = c("dataset", "condition"))

benchmarkingData_plot = benchmarkingData_plot %>% filter(dataset != "smartseq3_umis")

plot = ggplot(benchmarkingData_plot, aes(x = dataset, y = auroc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  #geom_jitter(data = benchmarkingData_plot, aes(x = dataset, y = auroc, group = condition, color = condition), position = position_dodge(width = 0.75)) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_smartseq3_AUROC_reads_0.08.svg'), width = 2.5, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_smartseq3_AUROC_reads_0.08.png'), width = 2.5, height = 3)


#### statistics

# reads (main figure), AUPRC

#benchmarkingData_formatted is ungrouped for this

benchmarkingData_reads = benchmarkingData_formatted %>% filter(dataset != "smartseq3_umis")

stats_auprc = benchmarkingData_reads %>%
  pivot_wider(names_from = condition, values_from = auprc, id_cols = c(dataset, sample, dbl_exp))

#DoubletFinder
stats_auprc_DoubletFinder = stats_auprc %>%
  select(dataset, DoubletFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = DoubletFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_DoubletFinder$TREX, stats_auprc_DoubletFinder$smartseq3_reads, paired = FALSE) # W = 157, p-value = 0.0388

#scDblFinder
stats_auprc_scDblFinder = stats_auprc %>%
  select(dataset, scDblFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = scDblFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_scDblFinder$TREX, stats_auprc_scDblFinder$smartseq3_reads, paired = FALSE) # W = 85, p-value = 0.3464

#hybrid
stats_auprc_hybrid = stats_auprc %>%
  select(dataset, hybrid) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = hybrid) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_hybrid$TREX, stats_auprc_hybrid$smartseq3_reads, paired = FALSE) # W = 9, p-value = 2.243e-06

#Scrublet
stats_auprc_Scrublet = stats_auprc %>%
  select(dataset, Scrublet) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = Scrublet) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_Scrublet$TREX, stats_auprc_Scrublet$smartseq3_reads, paired = FALSE) #W = 207, p-value = 2.243e-06


## AUROC


stats_auroc = benchmarkingData_reads %>%
  pivot_wider(names_from = condition, values_from = auroc, id_cols = c(dataset, sample, dbl_exp))


#DoubletFinder
stats_auroc_DoubletFinder = stats_auroc %>%
  select(dataset, DoubletFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = DoubletFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_DoubletFinder$TREX, stats_auroc_DoubletFinder$smartseq3_reads, paired = FALSE) # W = 144, p-value = 0.1304

#scDblFinder
stats_auroc_scDblFinder = stats_auroc %>%
  select(dataset, scDblFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = scDblFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_scDblFinder$TREX, stats_auroc_scDblFinder$smartseq3_reads, paired = FALSE) # W = 86, p-value = 0.3684

#hybrid
stats_auroc_hybrid = stats_auroc %>%
  select(dataset, hybrid) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = hybrid) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_hybrid$TREX, stats_auroc_hybrid$smartseq3_reads, paired = FALSE) # W = 87, p-value = 0.3913

#Scrublet
stats_auroc_Scrublet = stats_auroc %>%
  select(dataset, Scrublet) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = Scrublet) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_Scrublet$TREX, stats_auroc_Scrublet$smartseq3_reads, paired = FALSE) # W = 180, p-value = 0.001627



#### UMIs


benchmarkingData_umis = benchmarkingData_formatted %>% filter(dataset != "smartseq3_reads")

stats_auprc = benchmarkingData_umis %>%
  pivot_wider(names_from = condition, values_from = auprc, id_cols = c(dataset, sample, dbl_exp))

#DoubletFinder
stats_auprc_DoubletFinder = stats_auprc %>%
  select(dataset, DoubletFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = DoubletFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_DoubletFinder$TREX, stats_auprc_DoubletFinder$smartseq3_umis, paired = FALSE) # W = 205, p-value = 4.094e-05

#scDblFinder
stats_auprc_scDblFinder = stats_auprc %>%
  select(dataset, scDblFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = scDblFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_scDblFinder$TREX, stats_auprc_scDblFinder$smartseq3_umis, paired = FALSE) # W = 122, p-value = 0.5732

#hybrid
stats_auprc_hybrid = stats_auprc %>%
  select(dataset, hybrid) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = hybrid) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_hybrid$TREX, stats_auprc_hybrid$smartseq3_umis, paired = FALSE) # W = 36, p-value = 0.002283

#Scrublet
stats_auprc_Scrublet = stats_auprc %>%
  select(dataset, Scrublet) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = Scrublet) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_Scrublet$TREX, stats_auprc_Scrublet$smartseq3_umis, paired = FALSE) #W = 216, p-value = 2.312e-08


## AUROC


stats_auroc = benchmarkingData_umis %>%
  pivot_wider(names_from = condition, values_from = auroc, id_cols = c(dataset, sample, dbl_exp))


#DoubletFinder
stats_auroc_DoubletFinder = stats_auroc %>%
  select(dataset, DoubletFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = DoubletFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_DoubletFinder$TREX, stats_auroc_DoubletFinder$smartseq3_umis, paired = FALSE) # W = 205, p-value = 4.094e-05

#scDblFinder
stats_auroc_scDblFinder = stats_auroc %>%
  select(dataset, scDblFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = scDblFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_scDblFinder$TREX, stats_auroc_scDblFinder$smartseq3_umis, paired = FALSE) # W = 185, p-value = 0.0006576

#hybrid
stats_auroc_hybrid = stats_auroc %>%
  select(dataset, hybrid) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = hybrid) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_hybrid$TREX, stats_auroc_hybrid$smartseq3_umis, paired = FALSE) # W = 78, p-value = 0.2081

#Scrublet
stats_auroc_Scrublet = stats_auroc %>%
  select(dataset, Scrublet) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = dataset, values_from = Scrublet) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_Scrublet$TREX, stats_auroc_Scrublet$smartseq3_umis, paired = FALSE) #W = 216, p-value = 2.312e-08









  
