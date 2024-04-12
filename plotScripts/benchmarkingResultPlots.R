### plotting results of doublet  for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20240207
### Last edited by Madeline E Melzer on 20240319
# added updated paths for summation doublets and plotted, also re-plotted non-barcoded benchmarking data


library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(ggbreak)
library(scales)
library(glue)

set.seed = 23

plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/benchmarkingResults/"
dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/"

palette = c("classifier" = "#f59fff", 
            "scDblFinder" = "#b26c98", 
            "DoubletFinder" = "#feb325", 
            "Scrublet" = "#2474a9", 
            "hybrid" = "#03bfc4", 
            "scrambled" = "#007024",
            "otherMethods" = "#808285")

nonBarcodedDatasets = c("cline-ch", 
                "hm-12k", 
                "hm-6k", 
                "HMEC-orig-MULTI", 
                "HMEC-rep-MULTI",
                "HEK-HMEC-MULTI",
                #"J293t-dm", #the datas
                "mkidney-ch", 
                #"pbmc-1A-dm", 
                #"pbmc-1B-dm",
                #"pbmc-1C-dm",
                "pbmc-2ctrl-dm",
                "pbmc-2stim-dm",
                "pbmc-ch",
                "pdx-MULTI",
                "nuc-MULTI")

barcodedDatasets = c("FM01", 
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
                   "non" = "Jiang et al.", 
                   #"ClonMapper", 
                   "Biorxiv" = "Jain et al.", 
                   #"SPLINTR", 
                   #"TREX",
                   #"LARRY",
                   "cellTag" = "CellTag-multi",
                   "watermelon" = "Watermelon")

##########################################################################################################################
#                             
#  redownloaded the data to make sure I have finals and I am running it again here    
#                             
##########################################################################################################################

##########################################################################################################################
# final (FM01-FM08, non_cancer, SPLINTR, cellTag, TREX, watermelon, and smartseq3_umis)
##########################################################################################################################

dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/final/"

subfolders <- list.dirs(dataDirectory, recursive = FALSE)

final <- map_dfr(subfolders, ~{
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
#write.csv(final, paste0(dataDirectory, "allBenchmarking_fromFinalFolder.csv"))

final_datasets = final %>%
  filter(dataset %in% c("FM01", 
                        "FM02", 
                        "FM03", 
                        "FM04", 
                        "FM05", 
                        "FM06", 
                        "FM08", 
                        #"Biorxiv",
                        "non", #this is non_cancer, but naming was done improperly
                        #"ClonMapper", 
                        "SPLINTR", 
                        #"LARRY",
                        "cellTag",
                        "TREX",
                        "watermelon",
                        "smartseq3")) 
#need to filter it out for only the samples we want, 
filtered_final <- final_datasets %>%
  filter(
    !(
      (dataset == "SPLINTR" & sample %in% c("inVitro_KRAS", "retransplant")) |
        (dataset == "TREX" & sample == "brain1") |
        (dataset == "smartseq3" & sample %in% c("sample1", "reads_brain1", "reads_brain2"))
    )
  ) %>%
  mutate(dataset = str_replace(dataset, "non", "non_cancer")) %>%
  mutate(dataset = str_replace(dataset, "smartseq3", "smartseq3_umis"))
#write.csv(filtered_final, paste0(dataDirectory, "allBenchmarking_fromFinalFolder_filteredForUsableDatasets.csv"))



##########################################################################################################################
# non-barcoded
##########################################################################################################################

nonBarcoded <- final %>%
  filter (dataset %in% nonBarcodedDatasets)

#write.csv(nonBarcoded, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/final/nonBarcodedDatasets.csv")
# not used in final paper data- the data was all rerun. 

##########################################################################################################################
# partial_rerun (LARRY, ClonMapper, Biorxiv, smartseq3_reads)
##########################################################################################################################

subfolders <- list.dirs("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/partial_rerun/", recursive = FALSE) 

partial_rerun <- map_dfr(subfolders, ~{
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
partial_rerun = partial_rerun %>% mutate(dataset = str_replace(dataset, "smartseq3", "smartseq3_reads"))
#write.csv(partial_rerun, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/partial_rerun/allBenchmarking_fromPartialRerunFolder_allReranDatasetsUsable.csv")




##########################################################################################################################
# adding all datasets together
##########################################################################################################################

partial_rerun = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/partial_rerun/allBenchmarking_fromPartialRerunFolder_allReranDatasetsUsable.csv")
final = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/final/allBenchmarking_fromFinalFolder_filteredForUsableDatasets.csv")
#nonBarcoded = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/final/nonBarcodedDatasets.csv") # not used in final paper data
barcoded = bind_rows(final, partial_rerun)
#write_csv(barcoded, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/barcodedAllForPlot.csv")



barcodedDatasets = c("FM01", 
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


################################
#                              #
#    all datasets AUPRC        #
#                              #
################################

barcoded_formatted <- barcoded %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition)) %>%
  filter(dbl_act == 0.08) %>%
  group_by(dataset, condition, sample) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),
    dbl_act = mean(dbl_act, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp and dbl_act
    .groups = 'drop'      
  )

barcoded_summary = barcoded_formatted %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
barcoded_plot <- left_join(barcoded_formatted, barcoded_summary, by = c("dataset", "condition"))

barcoded_plot$dataset <- factor(barcoded_plot$dataset, levels = barcodedDatasets)
barcoded_plot <- barcoded_plot %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

# box whisker
plot = ggplot(barcoded_plot, aes(x = dataset, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUPRC_0.08.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUPRC_0.08.png'), width = 7, height = 4)

################################
#                              #
#    all datasets AUROC        #
#                              #
################################

barcoded_summary <- barcoded_formatted %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc)) # Replaced auprc with auroc

barcoded_plot <- left_join(barcoded_formatted, barcoded_summary, by = c("dataset", "condition"))

barcoded_plot$dataset <- factor(barcoded_plot$dataset, levels = barcodedDatasets)
barcoded_plot <- barcoded_plot %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

# Box whisker plot updated to use auroc instead of auprc
plot <- ggplot(barcoded_plot, aes(x = dataset, y = auroc, fill = condition)) + # Replaced auprc with auroc
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, color = condition), # Replaced auprc with auroc
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

# Saving the plot
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUROC_0.08.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUROC_0.08.png'), width = 7, height = 4)




################################
#                              #
#    barcoded/nonBarcoded      #
#                              #
################################



########## bringing together all of the updated nonBarcoded data from the rerun (20240320)
dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/withUpdatedNonBarcodedData"

subfolders <- list.dirs(dataDirectory, recursive = FALSE)

nonBarcoded_8 <- map_dfr(subfolders, ~{
  detection_file_path <- file.path(.x, "stats/all_detection_rates_original_benchmarking.tsv")
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
write.csv(nonBarcoded_8, glue("{dataDirectory}/nonBarcoded_8.csv"))


nonBarcoded_noSubsampling <- map_dfr(subfolders, ~{
  detection_file_path <- file.path(.x, "stats/all_detection_rates_original_benchmarking_no_mod.tsv")
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
write.csv(nonBarcoded_noSubsampling, glue("{dataDirectory}/nonBarcoded_noSubsampling.csv"))






### plotting

barcoded = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/barcodedAllForPlot.csv")
barcoded = barcoded %>% filter(dbl_act == 0.08)  %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition))
#write.csv(barcoded, glue("{dataDirectory}/barcoded_0.08_average.csv"))

barcoded_grouped = barcoded %>% group_by(dataset, condition, sample) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),
    dbl_act = mean(dbl_act, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp
    .groups = 'drop') %>% mutate(type = "average")
barcoded_grouped_filtered = barcoded_grouped %>% filter(dataset != "smartseq3_reads")  %>% filter(dataset != "smartseq3_umis")
write.csv(barcoded_grouped_filtered, glue("{dataDirectory}/barcoded_0.08_average_collapsedByDblExp.csv"))


nonBarcoded_8 = read.csv(glue("{dataDirectory}/nonBarcoded_8.csv"))
nonBarcoded_8 = nonBarcoded_8 %>% filter(dataset %in% nonBarcodedDatasets) %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition))
write.csv(nonBarcoded_8, glue("{dataDirectory}/nonBarcoded_0.08_filtered.csv"))



nonBarcoded_noSubsampling = read.csv(glue("{dataDirectory}/nonBarcoded_noSubsampling.csv"))
nonBarcoded_noSubsampling = nonBarcoded_noSubsampling %>% filter(dataset %in% nonBarcodedDatasets) %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition))
write.csv(nonBarcoded_noSubsampling, glue("{dataDirectory}/nonBarcoded_noSubsampling_filtered.csv"))




## for the 8% subsampled barcoded datasets

BCnBC = bind_rows(barcoded, nonBarcoded_8) %>% filter(dataset != "smartseq3_reads")  %>% filter(dataset != "smartseq3_umis")

BCnBC_formatted <- BCnBC %>%
  mutate(isBarcoded = case_when(
    dataset %in% barcodedDatasets    ~ "Barcoded",
    dataset %in% nonBarcodedDatasets ~ "NonBarcoded",
    TRUE                     ~ "Other"  # For datasets that don't fall into either category
  ))%>%
  dplyr::filter(isBarcoded %in% c("Barcoded", "NonBarcoded"))

### AUPRC

BCnBC_summary <- BCnBC_formatted %>%
  group_by(isBarcoded, condition) %>%
  summarise(
    sdAuprc = sd(auprc, na.rm = TRUE),  # Calculate standard deviation
    avgAuprc = mean(auprc, na.rm = TRUE),  # Calculate average
    .groups = 'drop'   # Drop grouping structure after summarizing
  ) 

BCnBC_plot = left_join(BCnBC_formatted, BCnBC_summary)

plot <- ggplot(BCnBC_plot, aes(x = isBarcoded, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = position_dodge(width = 0.9), width = 0.85) +
  geom_pointrange(data = BCnBC_summary, aes(x = isBarcoded, y = avgAuprc, ymin = avgAuprc - sdAuprc, ymax = avgAuprc + sdAuprc, color = condition, group = condition),
                  position = position_dodge(width = 0.9), shape = 16) +
  #geom_point(data = benchmarkingData_plot, aes(x = isBarcoded, y = auprc, group = condition, color = condition), position = position_dodge(width = 0.9)) +
  scale_color_manual(values = palette) +  # Ensure 'palette' is defined with your colors
  scale_fill_manual(values = palette) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUPRC_0.08.svg'), width = 2.5, height = 2.5)
ggsave(plot, file = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUPRC_0.08.png'), width = 2.5, height = 2.5)


### AUROC

BCnBC_summary <- BCnBC_formatted %>%
  group_by(isBarcoded, condition) %>%
  summarise(
    sdAuroc = sd(auroc, na.rm = TRUE),  # Calculate standard deviation
    avgAuroc = mean(auroc, na.rm = TRUE),  # Calculate average
    .groups = 'drop'   # Drop grouping structure after summarizing
  ) 

BCnBC_plot = left_join(BCnBC_formatted, BCnBC_summary, by = c("isBarcoded", "condition"))

plot <- ggplot(BCnBC_plot, aes(x = isBarcoded, y = auroc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = position_dodge(width = 0.9), width = 0.85) +
  geom_pointrange(data = BCnBC_summary, aes(x = isBarcoded, y = avgAuroc, ymin = avgAuroc - sdAuroc, ymax = avgAuroc + sdAuroc, color = condition, group = condition),
                  position = position_dodge(width = 0.9), shape = 16) +
  scale_color_manual(values = palette) +  # Ensure 'palette' is defined with your colors
  scale_fill_manual(values = palette) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, filename = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUROC_0.08.svg'), width = 2.5, height = 2.5)
ggsave(plot, filename = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUROC_0.08.png'), width = 2.5, height = 2.5)


## for the non subsampled barcoded datasets

BCnBC = bind_rows(barcoded, nonBarcoded_noSubsampling) %>% filter(dataset != "smartseq3_reads")  %>% filter(dataset != "smartseq3_umis")

BCnBC_formatted <- BCnBC %>%
  mutate(isBarcoded = case_when(
    dataset %in% barcodedDatasets    ~ "Barcoded",
    dataset %in% nonBarcodedDatasets ~ "NonBarcoded",
    TRUE                     ~ "Other"  # For datasets that don't fall into either category
  ))%>%
  dplyr::filter(isBarcoded %in% c("Barcoded", "NonBarcoded"))

### AUPRC

BCnBC_summary <- BCnBC_formatted %>%
  group_by(isBarcoded, condition) %>%
  summarise(
    sdAuprc = sd(auprc, na.rm = TRUE),  # Calculate standard deviation
    avgAuprc = mean(auprc, na.rm = TRUE),  # Calculate average
    .groups = 'drop'   # Drop grouping structure after summarizing
  ) 

BCnBC_plot = left_join(BCnBC_formatted, BCnBC_summary)

plot <- ggplot(BCnBC_plot, aes(x = isBarcoded, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = position_dodge(width = 0.9), width = 0.85) +
  geom_pointrange(data = BCnBC_summary, aes(x = isBarcoded, y = avgAuprc, ymin = avgAuprc - sdAuprc, ymax = avgAuprc + sdAuprc, color = condition, group = condition),
                  position = position_dodge(width = 0.9), shape = 16) +
  #geom_point(data = benchmarkingData_plot, aes(x = isBarcoded, y = auprc, group = condition, color = condition), position = position_dodge(width = 0.9)) +
  scale_color_manual(values = palette) +  # Ensure 'palette' is defined with your colors
  scale_fill_manual(values = palette) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUPRC_noSubsampling.svg'), width = 2.5, height = 2.5)
ggsave(plot, file = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUPRC_noSubsampling.png'), width = 2.5, height = 2.5)


### AUROC

BCnBC_summary <- BCnBC_formatted %>%
  group_by(isBarcoded, condition) %>%
  summarise(
    sdAuroc = sd(auroc, na.rm = TRUE),  # Calculate standard deviation
    avgAuroc = mean(auroc, na.rm = TRUE),  # Calculate average
    .groups = 'drop'   # Drop grouping structure after summarizing
  ) 

BCnBC_plot = left_join(BCnBC_formatted, BCnBC_summary, by = c("isBarcoded", "condition"))

plot <- ggplot(BCnBC_plot, aes(x = isBarcoded, y = auroc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = position_dodge(width = 0.9), width = 0.85) +
  geom_pointrange(data = BCnBC_summary, aes(x = isBarcoded, y = avgAuroc, ymin = avgAuroc - sdAuroc, ymax = avgAuroc + sdAuroc, color = condition, group = condition),
                  position = position_dodge(width = 0.9), shape = 16) +
  scale_color_manual(values = palette) +  # Ensure 'palette' is defined with your colors
  scale_fill_manual(values = palette) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, filename = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUROC_noSubsampling.svg'), width = 2.5, height = 2.5)
ggsave(plot, filename = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_AUROC_noSubsampling.png'), width = 2.5, height = 2.5)




### TNR for barcoded and non-barcoded 


tnrDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/TNR"

tnrData_barcoded = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/averagedAndSummedDoublets.csv")
tnrData_barcoded_filtered = tnrData_barcoded %>% filter(type == "average") %>% select(dataset, condition, sample, TNR) %>% mutate(isBarcoded = "Barcoded")

## formatting nonBarcoded TNR data:
tnrData_nonBarcoded = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/withUpdatedNonBarcodedData/TNR_8_benchmarking.csv") #note data directory

tnrData_nonBarcoded_filtered = tnrData_nonBarcoded %>% 
  dplyr::rename(condition = method, dataset = id) %>% 
  mutate(dataset = str_remove(dataset, "\\.rds$")) %>%
  filter(dataset %in% nonBarcodedDatasets) %>%
  mutate(isBarcoded = case_when(
    dataset %in% barcodedDatasets    ~ "Barcoded",
    dataset %in% nonBarcodedDatasets ~ "NonBarcoded",
    TRUE                     ~ "Other"  # For datasets that don't fall into either category
  )) %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition))


tnrData_BCnBC = bind_rows(tnrData_barcoded_filtered, tnrData_nonBarcoded_filtered) #length = 348 + 48 = 396
write.csv(tnrData_BCnBC, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/withUpdatedNonBarcodedData/TNR_8_barcodedNonBarcoded.csv") #note data directory


BCnBC_summary <- tnrData_BCnBC %>%
  group_by(isBarcoded, condition) %>%
  summarise(
    sdTNR = sd(TNR, na.rm = TRUE),  # Calculate standard deviation
    avgTNR = mean(TNR, na.rm = TRUE),  # Calculate average
    .groups = 'drop'   # Drop grouping structure after summarizing
  ) 

BCnBC_plot = left_join(tnrData_BCnBC, BCnBC_summary)

plot <- ggplot(BCnBC_plot, aes(x = isBarcoded, y = TNR, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = position_dodge(width = 0.9), width = 0.85) +
  geom_pointrange(data = BCnBC_summary, aes(x = isBarcoded, y = avgTNR, ymin = avgTNR - sdTNR, ymax = avgTNR + sdTNR, color = condition, group = condition),
                  position = position_dodge(width = 0.9), shape = 16) +
  #geom_point(data = benchmarkingData_plot, aes(x = isBarcoded, y = auprc, group = condition, color = condition), position = position_dodge(width = 0.9)) +
  scale_color_manual(values = palette) +  # Ensure 'palette' is defined with your colors
  scale_fill_manual(values = palette) +
  theme_classic() +
  scale_y_continuous(limits = c(0.8, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_TNR_0.08.svg'), width = 2.5, height = 2.5)
ggsave(plot, file = paste0(plotDirectory, 'barcodedNonBarcoded/benchmarkingResults_barcodedNonBarcoded_TNR_0.08.png'), width = 2.5, height = 2.5)

  



###### statistics

tnr = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/withUpdatedNonBarcodedData/TNR_8_barcodedNonBarcoded.csv")
nonBarcoded = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/withUpdatedNonBarcodedData/nonBarcoded_0.08_filtered.csv")
nonBarcoded = nonBarcoded %>% mutate(isBarcoded = "NonBarcoded")
barcoded = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/withUpdatedNonBarcodedData/barcoded_0.08_average_collapsedByDblExp.csv")
barcoded = barcoded %>% filter(dataset != "smartseq3_reads", dataset != "smartseq3_umis") %>% mutate(isBarcoded = "Barcoded")

barcodedNonBarcoded = bind_rows(barcoded, nonBarcoded) %>% 
  mutate(sample = str_replace(sample, "cancer_10kbarsiblingA", "10kbarsiblingA")) %>% 
  mutate(sample = str_replace(sample, "cancer_10kbarsiblingB", "10kbarsiblingB")) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename)) %>%
  select(dataset, sample, condition, auprc, auroc, dbl_act, isBarcoded) 

tnr_used = tnr %>% select(dataset, sample, condition, TNR)

allBarcodedNonBarcoded = inner_join(barcodedNonBarcoded, tnr_used, by = join_by(dataset, sample, condition))
write.csv(allBarcodedNonBarcoded, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv")


### averages:

barcodedForAvg = allBarcodedNonBarcoded %>% filter(isBarcoded == "Barcoded")

averages_by_condition <- barcodedForAvg %>%
  group_by(dataset) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),
    TNR = mean(TNR, na.rm = TRUE)
  ) %>%
  .groups = 'drop' %>%
  group_by(condition) %>%
  summarise(
    avg_auprc = mean(auprc, na.rm = TRUE),
    avg_auroc = mean(auroc, na.rm = TRUE),
    avg_TNR = mean(TNR, na.rm = TRUE)
  )






auprcData = allBarcodedNonBarcoded %>%
  pivot_wider(names_from = condition, values_from = auprc, id_cols = c(dataset, sample, isBarcoded))

# DoubletFinder
stats_auprc_DoubletFinder = auprcData %>%
  select(isBarcoded, DoubletFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = DoubletFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_DoubletFinder$Barcoded, stats_auprc_DoubletFinder$NonBarcoded, paired = FALSE) # W = 617, p-value = 0.311


# hybrid
stats_auprc_hybrid = auprcData %>%
  select(isBarcoded, hybrid) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = hybrid) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_hybrid$Barcoded, stats_auprc_hybrid$NonBarcoded, paired = FALSE) # W = 481, p-value = 0.6641


# scDblFinder
stats_auprc_scDblFinder = auprcData %>%
  select(isBarcoded, scDblFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = scDblFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_scDblFinder$Barcoded, stats_auprc_scDblFinder$NonBarcoded, paired = FALSE) # W = 51, p-value = 4.552e-07


# Scrublet
stats_auprc_Scrublet = auprcData %>%
  select(isBarcoded, Scrublet) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = Scrublet) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auprc_Scrublet$Barcoded, stats_auprc_Scrublet$NonBarcoded, paired = FALSE) # W = 571, p-value = 0.6031






######## AUROC Data


aurocData = allBarcodedNonBarcoded %>%
  pivot_wider(names_from = condition, values_from = auroc, id_cols = c(dataset, sample, isBarcoded))

# DoubletFinder
stats_auroc_DoubletFinder = aurocData %>%
  select(isBarcoded, DoubletFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = DoubletFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_DoubletFinder$Barcoded, stats_auroc_DoubletFinder$NonBarcoded, paired = FALSE) # W = 766, p-value = 0.009039


# hybrid
stats_auroc_hybrid = aurocData %>%
  select(isBarcoded, hybrid) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = hybrid) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_hybrid$Barcoded, stats_auroc_hybrid$NonBarcoded, paired = FALSE) # W = 733, p-value = 0.02402


# scDblFinder
stats_auroc_scDblFinder = aurocData %>%
  select(isBarcoded, scDblFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = scDblFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_scDblFinder$Barcoded, stats_auroc_scDblFinder$NonBarcoded, paired = FALSE) # W = 125, p-value = 2.129e-05


# Scrublet
stats_auroc_Scrublet = aurocData %>%
  select(isBarcoded, Scrublet) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = Scrublet) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_auroc_Scrublet$Barcoded, stats_auroc_Scrublet$NonBarcoded, paired = FALSE) # W = 749, p-value = 0.01517




### TNR
TNRData = allBarcodedNonBarcoded %>%
  pivot_wider(names_from = condition, values_from = TNR, id_cols = c(dataset, sample, isBarcoded))

# DoubletFinder
stats_TNR_DoubletFinder = TNRData %>%
  select(isBarcoded, DoubletFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = DoubletFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_TNR_DoubletFinder$Barcoded, stats_TNR_DoubletFinder$NonBarcoded, paired = FALSE) # W = 598, p-value = 0.4183

# hybrid
stats_TNR_hybrid = TNRData %>%
  select(isBarcoded, hybrid) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = hybrid) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_TNR_hybrid$Barcoded, stats_TNR_hybrid$NonBarcoded, paired = FALSE) # W = 373, p-value = 0.1114

# scDblFinder
stats_TNR_scDblFinder = TNRData %>%
  select(isBarcoded, scDblFinder) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = scDblFinder) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_TNR_scDblFinder$Barcoded, stats_TNR_scDblFinder$NonBarcoded, paired = FALSE) # W = 71, p-value = 1.366e-06

# Scrublet
stats_TNR_Scrublet = TNRData %>%
  select(isBarcoded, Scrublet) %>%
  mutate(row_id = row_number()) %>%                           # Add a unique row identifier
  pivot_wider(names_from = isBarcoded, values_from = Scrublet) %>%  # Pivot wider using the unique row identifier
  select(-row_id)
wilcox.test(stats_TNR_Scrublet$Barcoded, stats_TNR_Scrublet$NonBarcoded, paired = FALSE) # W = 643, p-value = 0.1964






















################################
#                              #
#    TNR                       #
#                              #
################################

tnrDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/TNR"

tnrData = read.csv(glue("{tnrDirectory}/TNR_8.csv"))

### formatting for plotting
tnrData <- tnrData %>% # Split 'id' column into multiple columns, ensure all parts are separated
  separate(id, into = c("dataset", "sample", "exp_act_part"), sep = "___", extra = "merge") %>%
  separate(exp_act_part, into = c("exp_part", "act_part"), sep = "__", extra = "merge")

tnrData <- tnrData %>% # Extract numeric values for 'dbl_exp' and 'dbl_act', adjusting the regular expression for 'dbl_act'
  mutate(dbl_exp = as.numeric(gsub("exp_([0-9.]+).*", "\\1", exp_part)),
         dbl_act = as.numeric(gsub("act_([0-9.]+).*", "\\1", act_part)))
tnrData <- tnrData %>% select(-exp_part, -act_part) # Drop the intermediate columns

#write.csv(tnrData, glue("{tnrDirectory}/TNR_plotted.csv"))
#tnrData = read.csv(glue("{tnrDirectory}/TNR_plotted.csv"))

### filtering so we only plot relevant samples

tnr_formatted = tnrData %>% 
  filter(dataset != "smartseq3")  %>% 
  filter(dataset != "smartseq3_reads")  %>% 
  filter(dataset != "smartseq3_umis") %>%
  filter(
    !(
      (dataset == "SPLINTR" & sample %in% c("inVitro_KRAS", "retransplant")) |
        (dataset == "TREX" & sample == "brain1") |
        (dataset == "smartseq3" & sample %in% c("sample1", "reads_brain1", "reads_brain2")) |
        (dataset == "LARRY" & sample %in% c("d4_nBC", "d4_R_4"))
    )
  ) %>%
  dplyr::rename(condition = method) %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition))

### plotting

## formatting TNR for average vs summation doublet comparison. This is average
tnr_collapsedByDblExp = tnr_formatted %>% group_by(dataset, condition, sample) %>%
  summarise(
    TNR = mean(TNR, na.rm = TRUE),
    dbl_act = mean(dbl_act, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp
    .groups = 'drop'
  ) %>% mutate(type = "average")

write.csv(tnr_collapsedByDblExp, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/TNR_barcoded_averageDoublets.csv")

tnr_summary <- tnr_formatted %>% group_by(dataset, condition) %>% summarise(sdTNR = sd(TNR), avgTNR = mean(TNR))

tnr_plot <- left_join(tnr_formatted, tnr_summary, by = c("dataset", "condition"))

tnr_plot$dataset <- factor(tnr_plot$dataset, levels = barcodedDatasets)
tnr_plot <- tnr_plot %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

plot <- ggplot(tnr_plot, aes(x = dataset, y = TNR, fill = condition)) + # Keep 'auroc' here if it's a column name in your data
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgTNR, ymin = avgTNR-sdTNR, ymax = avgTNR+sdTNR, color = condition),
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0.8,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, filename = paste0(plotDirectory, 'benchmarkingResults_barcoded_TNR_0.08.svg'), width = 7, height = 4)
ggsave(plot, filename = paste0(plotDirectory, 'benchmarkingResults_barcoded_TNR_0.08.png'), width = 7, height = 4)





















##########################################################################################################################
#                             
#  doublet detection methods on doublets formed via summation
#                             
##########################################################################################################################

dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/final_sum_rerun20240319/"

subfolders <- list.dirs(dataDirectory, recursive = FALSE)

final <- map_dfr(subfolders, ~{
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

dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/final_sum_rerun20240319/"
write.csv(final, paste0(dataDirectory, "allBenchmarking_sum.csv"))


#need to filter it out for only the samples we want, 
filtered_final <- final %>%
  filter(
    !(
      (dataset == "SPLINTR" & sample %in% c("inVitro_KRAS", "retransplant")) |
        (dataset == "TREX" & sample == "brain1") |
        (dataset == "smartseq3" & sample %in% c("sample1", "reads_brain1", "reads_brain2"))
    )
  ) %>%
  filter(dataset != "smartseq3") %>% #we dont need any of this data for summation analysis
  mutate(dataset = str_replace(dataset, "non", "non_cancer"))


write.csv(filtered_final, paste0(dataDirectory, "allBenchmarking_sum_filteredForUsableDatasets.csv"))


################################
#    all sum datasets AUPRC        
################################
dataDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/final_sum_rerun20240319/"
filtered_final = read.csv(paste0(dataDirectory, "allBenchmarking_sum_filteredForUsableDatasets.csv"))

barcoded = filtered_final #too lazy to rename from above (plotting the averaged doublets) so I am using this, here

barcoded_formatted <- barcoded %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition)) %>%
  filter(dbl_act == 0.08) %>%
  group_by(dataset, condition, sample) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),
    dbl_act = mean(dbl_act, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp and dbl_act
    .groups = 'drop'      
  )

#forSummedVsAveraged
barcoded_formatted_mutated = barcoded_formatted %>% mutate(type = "sum")
write.csv(barcoded_formatted_mutated, glue("{dataDirectory}barcoded_0.08_summed_collapsedByDblExp.csv"))





# auprc-specific
barcoded_summary = barcoded_formatted %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
barcoded_plot <- left_join(barcoded_formatted, barcoded_summary, by = c("dataset", "condition"))

barcoded_plot$dataset <- factor(barcoded_plot$dataset, levels = barcodedDatasets)
barcoded_plot <- barcoded_plot %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

# box whisker
plot = ggplot(barcoded_plot, aes(x = dataset, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUPRC_0.08_sum.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUPRC_0.08_sum.png'), width = 7, height = 4)


################################
#    all sum datasets AUROC        
################################

barcoded_summary <- barcoded_formatted %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc)) # Replaced auprc with auroc

barcoded_plot <- left_join(barcoded_formatted, barcoded_summary, by = c("dataset", "condition"))

barcoded_plot$dataset <- factor(barcoded_plot$dataset, levels = barcodedDatasets)
barcoded_plot <- barcoded_plot %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

# Box whisker plot updated to use auroc instead of auprc
plot <- ggplot(barcoded_plot, aes(x = dataset, y = auroc, fill = condition)) + # Replaced auprc with auroc
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, color = condition), # Replaced auprc with auroc
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

# Saving the plot
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUROC_0.08_sum.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcoded_AUROC_0.08_sum.png'), width = 7, height = 4)




################################
#    all sum datasets TNR        
################################

tnrDirectory = dataDirectory ="/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/final_sum_rerun20240319/"

tnrData = read.csv(glue("{tnrDirectory}/TNR_8_sum.csv"))

### formatting for plotting
tnrData <- tnrData %>% # Split 'id' column into multiple columns, ensure all parts are separated
  separate(id, into = c("dataset", "sample", "exp_act_part"), sep = "___", extra = "merge") %>%
  separate(exp_act_part, into = c("exp_part", "act_part"), sep = "__", extra = "merge")

tnrData <- tnrData %>% # Extract numeric values for 'dbl_exp' and 'dbl_act', adjusting the regular expression for 'dbl_act'
  mutate(dbl_exp = as.numeric(gsub("exp_([0-9.]+).*", "\\1", exp_part)),
         dbl_act = as.numeric(gsub("act_([0-9.]+).*", "\\1", act_part)))
tnrData <- tnrData %>% select(-exp_part, -act_part) # Drop the intermediate columns

### filtering so we only plot relevant samples

tnr_formatted = tnrData %>% 
  filter(dataset != "smartseq3")  %>% 
  filter(dataset != "smartseq3_reads")  %>% 
  filter(dataset != "smartseq3_umis") %>%
  filter(
    !(
      (dataset == "SPLINTR" & sample %in% c("inVitro_KRAS", "retransplant")) |
        (dataset == "TREX" & sample == "brain1") |
        (dataset == "smartseq3" & sample %in% c("sample1", "reads_brain1", "reads_brain2"))
    )
  ) %>%
  dplyr::rename(condition = method) %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition))

### plotting


## formatting TNR for average vs summation doublet comparison. This is summation
tnr_collapsedByDblExp = tnr_formatted %>% group_by(dataset, condition, sample) %>%
  summarise(
    TNR = mean(TNR, na.rm = TRUE),
    dbl_act = mean(dbl_act, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp
    .groups = 'drop'
  ) %>% mutate(type = "sum")
#write.csv(tnr_collapsedByDblExp, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/TNR_barcoded_summedDoublets.csv")


tnr_summary <- tnr_formatted %>% group_by(dataset, condition) %>% summarise(sdTNR = sd(TNR), avgTNR = mean(TNR))

tnr_plot <- left_join(tnr_formatted, tnr_summary, by = c("dataset", "condition"))

#write.csv(tnr_plot, glue("{tnrDirectory}/TNR_plotted_sum.csv"))

tnr_plot$dataset <- factor(tnr_plot$dataset, levels = barcodedDatasets)
tnr_plot <- tnr_plot %>%
  filter(dataset %in% barcodedDatasets) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

plot <- ggplot(tnr_plot, aes(x = dataset, y = TNR, fill = condition)) + # Keep 'auroc' here if it's a column name in your data
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgTNR, ymin = avgTNR-sdTNR, ymax = avgTNR+sdTNR, color = condition),
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0.8,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, filename = paste0(plotDirectory, 'benchmarkingResults_barcoded_TNR_0.08_sum.svg'), width = 7, height = 4)
ggsave(plot, filename = paste0(plotDirectory, 'benchmarkingResults_barcoded_TNR_0.08_sum.png'), width = 7, height = 4)























########################### to plot averaged vs summed, aggregate of all metrics (AUPRC, AUROC, TNR) 20240320

average = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/barcoded_0.08_average_collapsedByDblExp.csv")
sum = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/final_sum_rerun20240319/barcoded_0.08_summed_collapsedByDblExp.csv")

averageAndSum = bind_rows(average, sum) %>% 
  mutate(sample = str_replace(sample, "cancer_10kbarsiblingA", "10kbarsiblingA")) %>% 
  mutate(sample = str_replace(sample, "cancer_10kbarsiblingB", "10kbarsiblingB")) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename)) %>% 
  select(-X)

write.csv(averageAndSum,"/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/averagedAndSummedDoublets_AUPRC_AUROC.csv")

tnr_average = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/TNR_barcoded_averageDoublets.csv")
tnr_sum = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/TNR_barcoded_summedDoublets.csv")

averageAndSum_tnr = bind_rows(tnr_average, tnr_sum) %>% 
  select(-X)%>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

write.csv(averageAndSum_tnr,"/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/averagedAndSummedDoublets_TNR.csv")


averageAndSum_all = inner_join(averageAndSum, averageAndSum_tnr, by = c("dataset", "condition", "sample", "type", "dbl_act"))
write.csv(averageAndSum_all,"/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/averagedAndSummedDoublets.csv")



#### plotting for sanity
averageAndSum_all = averageAndSum_all %>% filter(type == "average")

barcoded_summary = averageAndSum_all %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
barcoded_plot <- left_join(averageAndSum_all, barcoded_summary, by = c("dataset", "condition"))

#barcoded_plot$dataset <- factor(barcoded_plot$dataset, levels = barcodedDatasets)
#barcoded_plot <- barcoded_plot# %>%
  #filter(dataset %in% barcodedDatasets)# %>%
  #mutate(dataset = recode(dataset, !!!barcodedRename))

# box whisker
plot = ggplot(barcoded_plot, aes(x = dataset, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot


























































































































################################################################################################################################################################################################################################################################################
# old
################################################################################################################################################################################################################################################################################


################################################################################################################################################################################################################################################################################
# old
################################################################################################################################################################################################################################################################################


################################################################################################################################################################################################################################################################################
# old
################################################################################################################################################################################################################################################################################


################################################################################################################################################################################################################################################################################
# old
################################################################################################################################################################################################################################################################################

dataset_sample_counts <- barcoded %>%
  group_by(dataset, sample) %>%
  summarise(row_count = n(), .groups = 'drop')

# View the counts
print(dataset_sample_counts)



allBenchmarking_final = bind_rows(data, filtered_data, doubletDetectionMethod_data_new)

write_csv(allBenchmarking_final, )


### all barcoded datasets
benchmarkingData_formatted <- benchmarkingData %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition),
         dataset = str_replace(dataset, "non", "non_cancer")) %>%
  group_by(dataset, condition, sample) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE), #average AUPRC for each sample, across all dbl_exp and dbl_act
    .groups = 'drop'      
  )









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
#write.csv(benchmarkingData, paste0(dataDirectory, "allBenchmarking.csv"))
#benchmarkingData = read_csv(paste0(dataDirectory, "allBenchmarking.csv"))


### all barcoded datasets
benchmarkingData_formatted <- benchmarkingData %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition),
         dataset = str_replace(dataset, "non", "non_cancer")) %>%
  group_by(dataset, condition, sample) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp and dbl_act
    .groups = 'drop'      
  )

#### filtering old benchmarking data, adding new benchmarking data

#for all except larry, clonmapper, biorxiv AND SMARTSEQ3, this is from "final benchmarking".
# The reason smartseq3 isn't here is because it had to be re-run. I did not use the re-run when I plotted for classifier, which is when I used "allBenchmarking.csv"
benchmarkingData = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/multiSampleCellIDResults/allBenchmarking.csv") %>% 
  select(-"...1") %>% 
  filter(dataset %in% c("FM01", 
                        "FM02", 
                        "FM03", 
                        "FM04", 
                        "FM05", 
                        "FM06", 
                        "FM08", 
                        #"Biorxiv",
                        "non", #this is non_cancer, but naming was done improperly
                        #"ClonMapper", 
                        "SPLINTR", 
                        #"LARRY",
                        "cellTag",
                        "TREX",
                        "watermelon")) 
#need to filter it out for only the samples we want
badSamples = c("SPLINTR" = "inVitro_KRAS", 
               "SPLINTR" = "retransplant", 
               "TREX" = "brain1")


subfolders <- list.dirs("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/LARRYClonMapperBiorxiv", recursive = FALSE) #this is with all of this morning's data. 

doubletDetectionMethod_data_new <- map_dfr(subfolders, ~{
  detection_file_path <- file.path(.x, "stats/all_detection_rates_0.1.tsv")
  if (file.exists(detection_file_path)) {
    read_tsv(detection_file_path) %>%
      filter(ID %in% names(id_mappings)) %>%
      mutate(dataset = id_mappings[ID],
             condition = basename(.x),
             auroc = roc,
             auprc = pr) %>%
      select(dataset, condition, auroc, auprc, dbl_act)
  } else {
    tibble()  # Return an empty tibble if the file doesn't exist
  }
})


#allData = bind_rows(data, filtered_data, doubletDetectionMethod_data_new)


allBenchmarking_final = bind_rows(data, filtered_data, doubletDetectionMethod_data_new)

write_csv(allBenchmarking_final, )


################################
#                              #
#    all datasets AUPRC        #
#                              #
################################

benchmarkingData_summary = benchmarkingData_formatted %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
benchmarkingData_plot <- left_join(benchmarkingData_formatted, benchmarkingData_summary, by = c("dataset", "condition"))

benchmarkingData_plot = benchmarkingData_plot %>% dplyr::filter(dataset %in% barcoded)
benchmarkingData_plot$dataset <- factor(benchmarkingData_plot$dataset, levels = barcoded)
benchmarkingData_plot <- benchmarkingData_plot %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

# box whisker
plot = ggplot(benchmarkingData_plot, aes(x = dataset, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUPRC.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUPRC.png'), width = 7, height = 4)

################################
#                              #
#    all datasets AUROC        #
#                              #
################################

benchmarkingData_summary = benchmarkingData_formatted %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc))
benchmarkingData_plot <- left_join(benchmarkingData_formatted, benchmarkingData_summary, by = c("dataset", "condition"))

benchmarkingData_plot = benchmarkingData_plot %>% dplyr::filter(dataset %in% barcoded)
benchmarkingData_plot$dataset <- factor(benchmarkingData_plot$dataset, levels = barcoded)
benchmarkingData_plot <- benchmarkingData_plot %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))


# box whisker
plot = ggplot(benchmarkingData_plot, aes(x = dataset, y = auroc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot


ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUROC.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUROC.png'), width = 7, height = 4)




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
#write.csv(benchmarkingData, paste0(dataDirectory, "allBenchmarking.csv"))
#benchmarkingData = read_csv(paste0(dataDirectory, "allBenchmarking.csv"))


### all barcoded datasets
benchmarkingData_formatted <- benchmarkingData %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition),
         dataset = str_replace(dataset, "non", "non_cancer")) %>%
  group_by(dataset, condition, sample) %>%
  summarise(
    auprc = mean(auprc, na.rm = TRUE),
    auroc = mean(auroc, na.rm = TRUE),#average AUPRC for each sample, across all dbl_exp and dbl_act
    .groups = 'drop'      
  )

#### filtering old benchmarking data, adding new benchmarking data

#for all except larry, clonmapper, biorxiv AND SMARTSEQ3, this is from "final benchmarking".
# The reason smartseq3 isn't here is because it had to be re-run. I did not use the re-run when I plotted for classifier, which is when I used "allBenchmarking.csv"
benchmarkingData = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/multiSampleCellIDResults/allBenchmarking.csv") %>% 
  select(-"...1") %>% 
  filter(dataset %in% c("FM01", 
                        "FM02", 
                        "FM03", 
                        "FM04", 
                        "FM05", 
                        "FM06", 
                        "FM08", 
                        #"Biorxiv",
                        "non", #this is non_cancer, but naming was done improperly
                        #"ClonMapper", 
                        "SPLINTR", 
                        #"LARRY",
                        "cellTag",
                        "TREX",
                        "watermelon")) 
#need to filter it out for only the samples we want
badSamples = c("SPLINTR" = "inVitro_KRAS", 
               "SPLINTR" = "retransplant", 
               "TREX" = "brain1")


subfolders <- list.dirs("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/LARRYClonMapperBiorxiv", recursive = FALSE) #this is with all of this morning's data. 

doubletDetectionMethod_data_new <- map_dfr(subfolders, ~{
  detection_file_path <- file.path(.x, "stats/all_detection_rates_0.1.tsv")
  if (file.exists(detection_file_path)) {
    read_tsv(detection_file_path) %>%
      filter(ID %in% names(id_mappings)) %>%
      mutate(dataset = id_mappings[ID],
             condition = basename(.x),
             auroc = roc,
             auprc = pr) %>%
      select(dataset, condition, auroc, auprc, dbl_act)
  } else {
    tibble()  # Return an empty tibble if the file doesn't exist
  }
})


#allData = bind_rows(data, filtered_data, doubletDetectionMethod_data_new)


allBenchmarking_final = bind_rows(data, filtered_data, doubletDetectionMethod_data_new)

write_csv(allBenchmarking_final, )


################################
#                              #
#    all datasets AUPRC        #
#                              #
################################

benchmarkingData_summary = benchmarkingData_formatted %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
benchmarkingData_plot <- left_join(benchmarkingData_formatted, benchmarkingData_summary, by = c("dataset", "condition"))

benchmarkingData_plot = benchmarkingData_plot %>% dplyr::filter(dataset %in% barcoded)
benchmarkingData_plot$dataset <- factor(benchmarkingData_plot$dataset, levels = barcoded)
benchmarkingData_plot <- benchmarkingData_plot %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

# box whisker
plot = ggplot(benchmarkingData_plot, aes(x = dataset, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUPRC.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUPRC.png'), width = 7, height = 4)

################################
#                              #
#    all datasets AUROC        #
#                              #
################################

benchmarkingData_summary = benchmarkingData_formatted %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc))
benchmarkingData_plot <- left_join(benchmarkingData_formatted, benchmarkingData_summary, by = c("dataset", "condition"))

benchmarkingData_plot = benchmarkingData_plot %>% dplyr::filter(dataset %in% barcoded)
benchmarkingData_plot$dataset <- factor(benchmarkingData_plot$dataset, levels = barcoded)
benchmarkingData_plot <- benchmarkingData_plot %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))


# box whisker
plot = ggplot(benchmarkingData_plot, aes(x = dataset, y = auroc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = 'dodge') +
  geom_pointrange(aes(x = dataset, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, color = condition), 
                  position = position_dodge(width = 0.75), size = 0.2, shape = 16) +
  scale_x_discrete() +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot


ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUROC.svg'), width = 7, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_allDatasets_AUROC.png'), width = 7, height = 4)


################################
#                              #
#    barcoded/nonBarcoded      #
#                              #
################################

benchmarkingData_formatted <- benchmarkingData %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"),
         condition = str_replace(condition, "scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         condition = as.factor(condition),
         dataset = str_replace(dataset, "non", "non_cancer")) %>%
  mutate(isBarcoded = case_when(
    dataset %in% barcoded    ~ "Barcoded",
    dataset %in% nonBarcoded ~ "NonBarcoded",
    TRUE                     ~ "Other"  # For datasets that don't fall into either category
  ))%>%
  dplyr::filter(isBarcoded %in% c("Barcoded", "NonBarcoded")) #%>%
#dplyr::filter(dbl_act == 0.08 | is.na(dbl_act))

benchmarkingData_summary <- benchmarkingData_formatted %>%
  group_by(isBarcoded, condition) %>%
  summarise(
    sdAuprc = sd(auprc, na.rm = TRUE),  # Calculate standard deviation
    avgAuprc = mean(auprc, na.rm = TRUE),  # Calculate average
    .groups = 'drop'   # Drop grouping structure after summarizing
  ) 

benchmarkingData_plot = left_join(benchmarkingData_formatted, benchmarkingData_summary)

plot <- ggplot(benchmarkingData_formatted, aes(x = isBarcoded, y = auprc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = position_dodge(width = 0.9), width = 0.85) +
  geom_pointrange(data = benchmarkingData_summary, aes(x = isBarcoded, y = avgAuprc, ymin = avgAuprc - sdAuprc, ymax = avgAuprc + sdAuprc, color = condition, group = condition),
                  position = position_dodge(width = 0.9), shape = 16) +
  #geom_point(data = benchmarkingData_plot, aes(x = isBarcoded, y = auprc, group = condition, color = condition), position = position_dodge(width = 0.9)) +
  scale_color_manual(values = palette) +  # Ensure 'palette' is defined with your colors
  scale_fill_manual(values = palette) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcodedNonBarcoded_AUPRC.svg'), width = 5, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcodedNonBarcoded_AUPRC.png'), width = 5, height = 4)


### AUROC

benchmarkingData_summary <- benchmarkingData_formatted %>%
  group_by(isBarcoded, condition) %>%
  summarise(
    sdAuroc = sd(auroc, na.rm = TRUE),  # Calculate standard deviation
    avgAuroc = mean(auroc, na.rm = TRUE),  # Calculate average
    .groups = 'drop'   # Drop grouping structure after summarizing
  ) 

benchmarkingData_plot = left_join(benchmarkingData_formatted, benchmarkingData_summary)

plot <- ggplot(benchmarkingData_formatted, aes(x = isBarcoded, y = auroc, fill = condition)) +
  geom_boxplot(aes(fill = condition), color = NA, alpha = 0.3, position = position_dodge(width = 0.9), width = 0.85) +
  geom_pointrange(data = benchmarkingData_summary, aes(x = isBarcoded, y = avgAuroc, ymin = avgAuroc - sdAuroc, ymax = avgAuroc + sdAuroc, color = condition, group = condition),
                  position = position_dodge(width = 0.9), shape = 16) +
  scale_color_manual(values = palette) +  # Ensure 'palette' is defined with your colors
  scale_fill_manual(values = palette) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcodedNonBarcoded_AUROC.svg'), width = 5, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'benchmarkingResults_barcodedNonBarcoded_AUROC.png'), width = 5, height = 4)

































