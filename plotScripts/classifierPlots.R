### plotting classifier figure for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20231220
### Last edited by Madeline E Melzer on 20240215

library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(ggbreak)
library(scales)

set.seed(23)

plot_directory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/classifier/"
dataset_dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/datasets/"

palette = c("classifier" = "#f59fff", 
            "scDblFinder" = "#b26c98", 
            "DoubletFinder" = "#feb325", 
            "Scrublet" = "#2474a9", 
            "hybrid" = "#03bfc4", 
            "scrambled" = "#007024",
            "otherMethods" = "#808285")


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
             #"s1nc_positiveControl")

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


# combine data
files <- list.files(dataset_dir, pattern = "\\.csv$", full.names = TRUE)
data <- files %>%
  map_dfr(~{
    df <- read_csv(.x)
    base_name <- tools::file_path_sans_ext(basename(.x))
    df_filtered <- filter(df, dataset == base_name)
    return(df_filtered)
  })

LARRYWatermelonCellTag = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/datasets/multiSampleCellIDsOLD/LARRYWatermelonCellTag.csv")
WatermelonCellTag = LARRYWatermelonCellTag %>% filter(dataset %in% c("watermelon", "cellTag"))

data = bind_rows(data, WatermelonCellTag)
#write.csv(data, paste0(dataset_dir, "combined/final/allDoubletRates_classifier.csv"))

# adding in results from other doublet detection
id_mappings <- c("Biorxiv_2_DMSO_B" = "Biorxiv",
                 "ClonMapper_FM1" = "ClonMapper",
                 "FM01_sample1" = "FM01",
                 "FM02_1-1uMPLX" = "FM02",
                 "FM03_DMSO-FM3-1uMPLX" = "FM03",
                 "FM04_BC18_B1" = "FM04",
                 "FM05_FM05-250nMPLX-A1" = "FM05",
                 "FM06_FM06-WM989Naive-1" = "FM06",
                 "FM08_run1_sample3" = "FM08",
                 "non_cancer_10kbarsiblingA" = "non_cancer",
                 "SPLINTR_chemoDay2_1" = "SPLINTR",
                 "TREX_brain2" = "TREX",
                 "LARRY_LK_d4_1" = "LARRY",
                 "watermelon_T47D_naive1",
                 "cellTag_d2-RNA-5" = "cellTag")

dataset_sample_pairs <- data.frame(
  dataset = c("Biorxiv", "ClonMapper", "FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "non", "SPLINTR", "TREX", "LARRY", "watermelon", "cellTag"),
  sample = c("1_DMSO_B", "FM1","sample1", "1-1uMPLX", "DMSO-FM3-1uMPLX", "BC18_B1", "FM05-250nMPLX-A1", "FM06-WM989Naive-1", "run1_sample3", "cancer_10kbarsiblingA", "chemoDay2_1", "brain2", "LK_d4_1", "T47D-naive-1", "d2-RNA-5")
)

# subfolders <- list.dirs("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/multiSampleCellIDResults", recursive = FALSE)
# 
# doubletDetectionMethod_data <- map_dfr(subfolders, ~{
#   detection_file_path <- file.path(.x, "all_detection_rates.tsv")
#   
#   if (file.exists(detection_file_path)) {
#     read_tsv(detection_file_path) %>%
#       filter(ID %in% names(id_mappings)) %>%
#       mutate(dataset = id_mappings[ID],
#              condition = basename(.x),
#              auroc = roc,
#              auprc = pr) %>%
#       select(dataset, condition, auroc, auprc, dbl_act)
#   } else {
#     tibble()  # Return an empty tibble if the file doesn't exist
#   }
# })


#for all except larry, clonmapper, biorxiv
doubletDetectionMethod_data = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/multiSampleCellIDResults/allBenchmarking.csv") %>% 
  select(-"...1") %>% 
  filter(dbl_act == 0.10) %>%
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
filtered_data <- doubletDetectionMethod_data %>%
  inner_join(dataset_sample_pairs, by = c("dataset", "sample"))



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


allData = bind_rows(data, filtered_data, doubletDetectionMethod_data_new)


#data = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/datasets/combined/allDatasets_classifier.csv")
#data2 = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/datasets/LARRYWatermelonCellTag.csv")
#data_combined <- data %>% bind_rows(doubletDetectionMethod_data) %>% filter(dbl_act == 0.10)
#data_combined = data_combined %>% bind_rows(data2)


#write.csv(allData, paste0(dataset_dir, "combined/final/allClassifierData.csv")) #20240215, 06:25 pm
allData = read_csv(paste0(dataset_dir, "combined/final/allClassifierData.csv"))
#allData = allData %>% mutate(dataset = str_replace(dataset, "^non$", "non_cancer")) #did once so do not have to do again

#data = read_csv("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/datasets/allDatasetsAndOtherMethods.csv")
data_formatted <- allData %>%
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "classifier: doublet control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "classifier: singlet control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
         #condition = str_replace(condition, "s1nc_positiveControl", "combined data"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled gene position"),
         condition = if_else(condition == as.character(dataset), "classifier", condition),
         condition = str_replace(condition, "0.10", "classifier"),
         condition = str_replace(condition, "0.20", "classifier"),
         condition = str_replace(condition, "FM01_params", "FM01 optimized parameters"),
         condition = str_replace(condition, "DoubletFinder", "DoubletFinder"),
         condition = str_replace(condition, "Scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"))

data_plot = data_formatted %>% 
  filter(dataset %in% barcoded) %>%
  #filter(condition %in% c("classifier", "scrambled", "DoubletFinder", "scDblFinder", "Scrublet", "hybrid"))  %>%
  mutate(dataset = factor(dataset, levels = barcoded)) %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

data_plot_toSave = data_plot %>% select(-"...1", -"...2", -auroc.groups)

#write_csv(data_plot_toSave, paste0(dataset_dir, "combined/final/allClassifierData_formatted.csv")) #20240325

## statistics 20240321

averages_by_condition <- data_plot %>%
  group_by(condition) %>%
  summarise(
    avg_auprc = mean(auprc, na.rm = TRUE),
    avg_auroc = mean(auroc, na.rm = TRUE)
  )



classifierStats = data_plot %>%
  mutate(condition = str_replace(condition, "DoubletFinder", "otherMethods"),
         condition = str_replace(condition, "scDblFinder", "otherMethods"),
         condition = str_replace(condition, "Scrublet", "otherMethods"),
         condition = str_replace(condition, "hybrid", "otherMethods"))

auprc_classifier <- classifierStats %>%
  filter(condition == "classifier") %>%
  pull(auprc)
auprc_otherMethods <- classifierStats %>%
  filter(condition == "otherMethods") %>%
  pull(auprc)
auprc_scrambled <- classifierStats %>%
  filter(condition == "scrambled") %>%
  pull(auprc)

# Perform the Wilcoxon rank-sum test
wilcox.test(auprc_classifier, auprc_otherMethods, paired = FALSE) #W = 53982, p-value < 2.2e-16
wilcox.test(auprc_classifier, auprc_scrambled, paired = FALSE) #W = 22500, p-value < 2.2e-16
t.test(auprc_classifier, auprc_otherMethods, paired = FALSE) #t = 71.036, df = 488.88, p-value < 2.2e-16
t.test(auprc_classifier, auprc_scrambled, paired = FALSE) #t = 191.23, df = 212.95, p-value < 2.2e-16

auroc_classifier <- classifierStats %>%
  filter(condition == "classifier") %>%
  pull(auroc)
auroc_otherMethods <- classifierStats %>%
  filter(condition == "otherMethods") %>%
  pull(auroc)
auroc_scrambled <- classifierStats %>%
  filter(condition == "scrambled") %>%
  pull(auroc)

# Perform the Wilcoxon rank-sum test
wilcox.test(auroc_classifier, auroc_otherMethods, paired = FALSE) #W = 53950, p-value < 2.2e-16
wilcox.test(auroc_classifier, auroc_scrambled, paired = FALSE) #W = 53950, p-value < 2.2e-16
t.test(auroc_classifier, auroc_otherMethods, paired = FALSE) #t = 31.63, df = 377.38, p-value < 2.2e-16
wilcox.test(auroc_classifier, auroc_scrambled, paired = FALSE) #t = 139.64, df = 187.48, p-value < 2.2e-16


publicData_classifier = read.csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/Revisions/plotData/classifier/allClassifierData_formatted.csv")
data_plot = publicData_classifier %>% select(-sdAuprc, -avgAuprc) %>% mutate(dataset = factor(dataset, levels = barcodedDatasets_fromPublicData))



### AUPRC
data_plot_summary = data_plot %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc), auroc.groups = 'drop')
data_plot <- left_join(data_plot, data_plot_summary, by = c("dataset", "condition"))
plot = ggplot(data_plot, aes(x = dataset, y = auprc)) +
  geom_ribbon(data = data_plot, aes(ymin = avgAuprc - sdAuprc, ymax = avgAuprc + sdAuprc, group = condition, fill = condition), alpha = 0.5) +
  #geom_ribbon(data = filter(data_plot, condition == "hybrid"),
  #            aes(ymin = avgAuprc - 0.002, ymax = avgAuprc + 0.002, fill = condition, group = condition), alpha = 0.5) +
  #geom_jitter(data = data_plot, aes(x=dataset, y=auprc),size=1, shape = 16, width = 0.3, height = 0.01) +
  #geom_pointrange(data = data_plot, aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition)) +
  geom_point(data = data_plot, aes(x = dataset, y =avgAuprc, color = condition), size = 1.5, shape = 16) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1.04))+ #, breaks = c(0.25, 0.50, 0.75, 1.0)) +  # Set custom y-axis breaks
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_allDatasets.svg'), width = 4.5, height = 3)
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_allDatasets.png'), width = 4.5, height = 3)

### AUROC
data_plot_summary = data_plot %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc), auroc.groups = 'drop')
data_plot <- left_join(data_plot, data_plot_summary, by = c("dataset", "condition"))
plot = ggplot(data_plot, aes(x = dataset, y = auroc)) +
  geom_ribbon(data = data_plot, aes(ymin = avgAuroc - sdAuroc, ymax = avgAuroc + sdAuroc, group = condition, fill = condition), alpha = 0.5) +
  #geom_ribbon(data = filter(data_plot, condition == "hybrid"),
  #            aes(ymin = avgAuprc - 0.002, ymax = avgAuprc + 0.002, fill = condition, group = condition), alpha = 0.5) +
  #geom_jitter(data = data_plot, aes(x=dataset, y=auprc),size=1, shape = 16, width = 0.3, height = 0.01) +
  #geom_pointrange(data = data_plot, aes(x = dataset, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, color = condition)) +
  geom_point(data = data_plot, aes(x = dataset, y =avgAuroc, color = condition), size = 1.5, shape = 16) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1.04))+ #, breaks = c(0.25, 0.50, 0.75, 1.0)) +  # Set custom y-axis breaks
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plot_directory, 'AUROC_allDatasets.svg'), width = 4.5, height = 3)
#ggsave(plot, file = paste0(plot_directory, 'AUROC_allDatasets.png'), width = 4.5, height = 3)

################################## allDatasets, in brief
publicData_classifier = read.csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/Revisions/plotData/classifier/allClassifierData_formatted.csv")
allData = publicData_classifier %>% select(-sdAuprc, -avgAuprc) %>% mutate(dataset = factor(dataset, levels = barcodedDatasets_fromPublicData))

data_formatted <- allData %>%
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "doublet control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "singlet control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled"),
         condition = str_replace(condition, "DoubletFinder", "otherMethods"),
         condition = str_replace(condition, "Scrublet", "otherMethods"),
         condition = str_replace(condition, "scDblFinder", "otherMethods"),
         condition = str_replace(condition, "hybrid", "otherMethods"),
         condition = if_else(condition == as.character(dataset), "classifier", condition))

data_plot = data_formatted %>% 
  #filter(dataset %in% barcoded) %>%
  filter(condition %in% c("classifier", "scrambled", "otherMethods"))
data_plot_summary = data_plot %>% group_by(condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
data_plot <- left_join(data_plot, data_plot_summary, by = "condition")
data_plot$condition <- fct_relevel(data_plot$condition, "scrambled", "otherMethods", "classifier")
plot = ggplot(data_plot, aes(x = condition, y = auprc)) +
  #geom_violin(data = data_plot, aes(x = condition, y = auprc, fill = condition)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auprc, color = condition), size=1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = condition, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), color = "black", size = 1) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) + # making sure the error bar fits
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_allDatasets_sum.svg'), width = 3, height = 6)
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_allDatasets_sum.png'), width = 3, height = 6)


data_plot = data_formatted %>% 
  #filter(dataset %in% barcoded) %>%
  filter(condition %in% c("classifier", "scrambled", "otherMethods"))
data_plot_summary = data_plot %>% group_by(condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc))
data_plot <- left_join(data_plot, data_plot_summary, by = "condition")
data_plot$condition <- fct_relevel(data_plot$condition, "scrambled", "otherMethods", "classifier")
plot = ggplot(data_plot, aes(x = condition, y = auroc)) +
  #geom_violin(data = data_plot, aes (x = condition, y = auprc, fill = condition)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auroc, color = condition),size= 1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = condition, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, fill = "condition"), color = "black", size = 1) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) + #making sure the error bar fits
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plot_directory, 'AUROC_allDatasets_sum.svg'), width = 3, height = 6)
#ggsave(plot, file = paste0(plot_directory, 'AUROC_allDatasets_sum.png'), width = 3, height = 6)




auprc_classifier <- data_plot %>%
  filter(condition == "classifier") %>%
  pull(auprc)
auprc_otherMethods <- data_plot %>%
  filter(condition == "otherMethods") %>%
  pull(auprc)

# Perform the Wilcoxon rank-sum test
p_value <- wilcox.test(auprc_classifier, auprc_otherMethods)$p.value
print(p_value)



################################## variable doublet rate

dataset_dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/variableDoubletRates"

# adding in results from other doublet detection methods
files <- list.files(dataset_dir, pattern = "\\.csv$", full.names = TRUE)
data <- files %>%
  map_dfr(~{
    df <- read_csv(.x)
    base_name <- tools::file_path_sans_ext(basename(.x))
    base_name = as.numeric(base_name)
    df_filtered <- filter(df, dataset == base_name)
    return(df_filtered)
  })

#write_csv(data, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/variableDoubletRates/allClassifierVariableData.csv")
publicData_variableClassifier = read.csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/Revisions/plotData/classifier/allClassifierVariableData.csv")
data = publicData_variableClassifier

id_mappings <- c("FM01_sample1" = "FM01")

subfolders <- list.dirs("/Users/mem3579/quest/ZhangMelzerEtAl/data/detectionMethods/", recursive = FALSE)

doubletDetectionMethod_data = read_csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/Revisions/plotData/benchmarking/barcodedAllForPlot.csv") %>% 
  #select(-"...1") %>% 
  filter(dataset %in% c("FM01")) %>%
  filter(sample %in% c("sample1"))


doubletDetectionMethod_data_mutated = doubletDetectionMethod_data %>% mutate(dataset = dbl_act) %>% 
  mutate(dataset = as.character(dataset)) %>%
  select(dataset, condition, auprc, auroc)

# doubletDetectionMethod_data <- map_dfr(subfolders, ~{
#   detectionMethod = basename(.x)
#   detection_file_paths <- list.files(.x, pattern = "all_detection_rates_0.*\\.tsv", full.names = TRUE)
#   
#   map_dfr(detection_file_paths, ~{
#     if (file.exists(.x)) {
#       read_tsv(.x) %>%
#         filter(ID %in% names(id_mappings)) %>%
#         mutate(dataset = str_extract(basename(.x), "0\\.\\d+"),  # Extract detection rate like "0.1", "0.2", etc.
#                condition = detectionMethod,
#                auroc = roc,
#                auprc = pr) %>%
#         select(dataset, condition, auroc, auprc)
#     } else {
#       tibble()  # Return an empty tibble if the file doesn't exist
#     }
#   })
# })

data <- data %>% mutate(dataset = as.character(dataset)) %>%
  select(dataset, condition, auprc, auroc)
data_combined <- data %>% bind_rows(doubletDetectionMethod_data_mutated)
#write_csv(data_combined, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/variableDoubletRates/allVariableData.csv")
data = data_combined

data_formatted <- data %>%
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "classifier: doublet control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "classifier: singlet control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
         dataset = str_replace(dataset, "\\b0.1\\b", "0.10"),
         dataset = str_replace(dataset, "\\b0.2\\b", "0.20"),
         #condition = str_replace(condition, "s1nc_positiveControl", "combined data"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled gene position"),
         condition = if_else(condition == as.character(dataset), "classifier", condition),
         condition = str_replace(condition, "0.10", "classifier"),
         condition = str_replace(condition, "0.20", "classifier"),
         condition = str_replace(condition, "FM01_params", "FM01 optimized parameters"),
         condition = str_replace(condition, "DoubletFinder", "DoubletFinder"),
         condition = str_replace(condition, "Scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"))

### AUPRC
data_plot = data_formatted %>% 
  filter(dataset %in% c("0.05", "0.08", "0.10", "0.15", "0.20", "0.25")) %>%
  filter(condition %in% c("classifier", "scrambled", "DoubletFinder", "scDblFinder", "Scrublet", "hybrid"))
#data_plot$dataset <- fct_relevel(data_plot$dataset, "s1nc_positiveControl", "FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "non_cancer", "Biorxiv", "TREX", "SPLINTR", "ClonMapper")
data_plot_summary = data_plot %>% group_by(dataset, condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
data_plot <- left_join(data_plot, data_plot_summary, by = c("dataset", "condition"))
plot <- ggplot(data_plot, aes(x = dataset, y = avgAuprc, group = condition, color = condition)) +
  geom_line() +
  geom_ribbon(aes(ymin = avgAuprc - sdAuprc, ymax = avgAuprc + sdAuprc, fill = condition), alpha = 0.5, colour = NA) +  # Set ribbon border to transparent
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_variableDoubletRates.svg'), width = 5, height = 4)
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_variableDoubletRates.png'), width = 5, height = 4)


### AUROC

data_plot = data_formatted %>% 
  filter(dataset %in% c("0.05", "0.08", "0.10", "0.15", "0.20", "0.25")) %>%
  filter(condition %in% c("classifier", "scrambled", "DoubletFinder", "scDblFinder", "Scrublet", "hybrid"))
#data_plot$dataset <- fct_relevel(data_plot$dataset, "s1nc_positiveControl", "FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "non_cancer", "Biorxiv", "TREX", "SPLINTR", "ClonMapper")
data_plot_summary = data_plot %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc))
data_plot <- left_join(data_plot, data_plot_summary, by = c("dataset", "condition"))
plot <- ggplot(data_plot, aes(x = dataset, y = avgAuroc, group = condition, color = condition)) +
  geom_line() +
  geom_ribbon(aes(ymin = avgAuroc - sdAuroc, ymax = avgAuroc + sdAuroc, fill = condition), alpha = 0.5, colour = NA) +  # Set ribbon border to transparent
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot
#ggsave(plot, file = paste0(plot_directory, 'AUROC_variableDoubletRates.svg'), width = 5, height = 4)
#ggsave(plot, file = paste0(plot_directory, 'AUROC_variableDoubletRates.png'), width = 5, height = 4)



############################### sample 1 sample 2 AUPRC

dataset_dir = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/datasets/"

#note: for these files, the dataset is the dataset being classified, the condition is the dataset which the classifier was trained on.
#ex. dataset sample2 and condition sample1 means the classifier trained on sample 1 is being used to classify sample2
s1d = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/summary_s1_bothClassifiers.csv")
s2d = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/summary_s2_bothClassifiers.csv")
s1s2d = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/summary_s1s2_bothClassifiers.csv")

data = s1d %>% bind_rows(s2d)

# adding in results from other doublet detection methods

id_mappings <- c("FM01_sample1" = "sample1",
                 "FM01_sample2" = "sample2")

subfolders <- list.dirs("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/multiSampleCellIDResults", recursive = FALSE)

doubletDetectionMethod_data <- map_dfr(subfolders, ~{
  detection_file_path <- file.path(.x, "stats/all_detection_rates_0.1.tsv")
  
  if (file.exists(detection_file_path)) {
    read_tsv(detection_file_path) %>%
      filter(ID %in% names(id_mappings)) %>%
      mutate(dataset = id_mappings[ID],
             condition = basename(.x),
             auroc = roc,
             auprc = pr) %>%
      select(dataset, condition, auroc, auprc)
  } else {
    tibble()  # Return an empty tibble if the file doesn't exist
  }
})



data_combined <- data %>% bind_rows(doubletDetectionMethod_data)
data = data_combined
#write_csv(data, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/alls1s2ForPlot.csv")
data = read.csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/Revisions/plotData/classifier/alls1s2ForPlot_formatted.csv")

# formatting and plotting
data_formatted <- data %>%
  filter(dataset != condition) %>% #eliminating the self-classifier 
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "doublet control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "singlet control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled"),
         condition = str_replace(condition, "DoubletFinder", "DoubletFinder"),
         condition = str_replace(condition, "Scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"))
         #condition = str_replace(condition, "sample1", "classifier"),
         #condition = str_replace(condition, "sample2", "classifier"))


#write_csv(data_formatted, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/alls1s2ForPlot_formatted.csv") #20240325

data_summary = data_formatted %>% group_by(condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc), auprc.groups = 'drop')
data_joined <- left_join(data_formatted, data_summary, by = "condition")
data_plot = data_joined %>% filter(condition %in% c("classifier", "DoubletFinder", "scDblFinder", "Scrublet", "hybrid", "scrambled"))
plot = ggplot(data_plot, aes(x = dataset, y = auprc)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auprc, color = condition),size=1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = condition, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), color = "black", size = 1) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_sample1Sample2ClassifiersAndMethods.svg'), width = 6, height = 4)
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_sample1Sample2ClassifiersAndMethods.png'), width = 6, height = 4)

### with algorithms combined
data_formatted <- data %>%
  filter(dataset != condition) %>% #eliminating the self-classifier 
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "doublet control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "singlet control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled"),
         condition = str_replace(condition, "DoubletFinder", "otherMethods"),
         condition = str_replace(condition, "Scrublet", "otherMethods"),
         condition = str_replace(condition, "scDblFinder", "otherMethods"),
         condition = str_replace(condition, "hybrid", "otherMethods"),
         condition = str_replace(condition, "sample1", "classifier"),
         condition = str_replace(condition, "sample2", "classifier"))

### AUPRC
data_summary = data_formatted %>% group_by(condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc), auprc.groups = 'drop')
data_joined <- left_join(data_formatted, data_summary, by = "condition")
data_plot = data_joined %>% filter(condition %in% c("otherMethods", "classifier"))
data_plot$condition <- fct_relevel(data_plot$condition, "otherMethods", "classifier")
plot = ggplot(data_plot, aes(x = dataset, y = auprc)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auprc, color = condition), size=1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = condition, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), color = "black", size = 1) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_sample1Sample2ClassifiersAndMethods_sum.svg'), width = 3, height = 6)
#ggsave(plot, file = paste0(plot_directory, 'AUPRC_sample1Sample2ClassifiersAndMethods_sum.png'), width = 3, height = 6)


### AUROC
data_summary = data_formatted %>% group_by(condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc), auroc.groups = 'drop')
data_joined <- left_join(data_formatted, data_summary, by = "condition")
data_plot = data_joined %>% filter(condition %in% c("otherMethods", "classifier"))
data_plot$condition <- fct_relevel(data_plot$condition, "otherMethods", "classifier")
plot = ggplot(data_plot, aes(x = dataset, y = auroc)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auroc, color = condition), size=1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = condition, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, fill = "condition"), color = "black", size = 1) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot
#ggsave(plot, file = paste0(plot_directory, 'AUROC_sample1Sample2ClassifiersAndMethods_sum.svg'), width = 3, height = 6)
#ggsave(plot, file = paste0(plot_directory, 'AUROC_sample1Sample2ClassifiersAndMethods_sum.png'), width = 3, height = 6)





# Extract 'auprc' values for each condition into separate vectors
auprc_classifier <- data_plot %>% #mean 0.6873221
  filter(condition == "classifier") %>%
  pull(auprc)
auprc_otherMethods <- data_plot %>% #mean 0.3571271
  filter(condition == "otherMethods") %>%
  pull(auprc)

# Perform the Wilcoxon rank-sum test
wilcox.test(auprc_classifier, auprc_otherMethods, paired = FALSE) #W = 9482, p-value < 2.2e-16

auroc_classifier <- data_plot %>% #mean 0.9180228
  filter(condition == "classifier") %>%
  pull(auroc)
auroc_otherMethods <- data_plot %>% #mean 0.779673
  filter(condition == "otherMethods") %>%
  pull(auroc)

# Perform the Wilcoxon rank-sum test
wilcox.test(auroc_classifier, auroc_otherMethods, paired = FALSE) #W = 9562, p-value < 2.2e-16







### for plotting more granular classifier data

s1s2d = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/summary_s1s2_bothClassifiers.csv")
data_withBoth = data %>% bind_rows(s1s2d)
#write.csv(data_withBoth, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/summary_s1s2_bothClassifiers_andSelf.csv")
data_withBoth = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/results/datasetSummaries/s1s2/summary_s1s2_bothClassifiers_andSelf.csv")
data_formatted <- data_withBoth %>%
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "doublet control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "singlet control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled"),
         condition = str_replace(condition, "DoubletFinder", "DoubletFinder"),
         condition = str_replace(condition, "Scrublet", "Scrublet"),
         condition = str_replace(condition, "scDblFinder", "scDblFinder"),
         dataset = str_replace(dataset, "s1s2", "sample1 and sample2")) %>%
  select(-"best_params")


data_summary = data_formatted %>% group_by(dataset, condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc))
data_joined <- left_join(data_formatted, data_summary, by = c("dataset", "condition"))
data_plot = data_joined %>% filter(condition %in% c("sample1", "sample2"))
sample_palette = c("sample1" = "#00888d", "sample2" = "#ff87a9")
plot = ggplot(data_plot, aes(x = dataset, y = auroc)) +
  geom_jitter(data = data_plot, aes(x=dataset, y=auroc, color = condition),size=1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = dataset, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, fill = condition), size = 1) +
  scale_color_manual(values = sample_palette) + 
  scale_fill_manual(values = sample_palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0.25, 0.50, 0.75, 1.0)) +  # Set custom y-axis breaks
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot
#ggsave(plot, file = paste0(plot_directory, 'AUROC_sample1Sample2CombinedClassifiers_sum.svg'), width = 3, height = 6)
#ggsave(plot, file = paste0(plot_directory, 'AUROC_sample1Sample2CombinedClassifiers_sum.png'), width = 3, height = 6)


################################
levels(T47D_seurat$active.ident)

##########################################################################################################################
# supplements with controls
##########################################################################################################################

#allData = read_csv(paste0(dataset_dir, "combined/final/allClassifierData.csv")) #20240215, 1:01 am

data_formatted <- allData %>%
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "doublet only control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "singlet only control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled control"),
         condition = str_replace(condition, "s1nc_positiveControl", "positive control"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled control"),
         condition = if_else(condition == as.character(dataset), "classifier", condition),
         condition = str_replace(condition, "FM01_params", "FM01 optimized parameters"),
         condition = str_replace(condition, "DoubletFinder", "otherMethods"),
         condition = str_replace(condition, "Scrublet", "otherMethods"),
         condition = str_replace(condition, "scDblFinder", "otherMethods"),
         condition = str_replace(condition, "hybrid", "otherMethods"))

### AUPRC
data_plot = data_formatted %>% 
  filter(dataset %in% barcoded) %>%
  filter(condition %in% c("classifier", "doublet only control", "singlet only control", "positive control", "shuffled control", "scrambled control", "otherMethods"))
data_plot_summary = data_plot %>% group_by(condition) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc))
data_plot <- left_join(data_plot, data_plot_summary, by = "condition")
data_plot$condition <- fct_relevel(data_plot$condition,  "scrambled control", "singlet only control", "doublet only control" ,"otherMethods", "shuffled control", "positive control", "classifier")

plot = ggplot(data_plot, aes(x = condition, y = auprc)) +
  #geom_violin(data = data_plot, aes(x = condition, y = auprc, fill = condition)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auprc, color = condition), size=1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = condition, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), color = "black", size = 1) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) + # making sure the error bar fits
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plot_directory, 'classifierControls_AUPRC.svg'), width = 8, height = 6)
#ggsave(plot, file = paste0(plot_directory, 'classifierControls_AUPRC.png'), width = 8, height = 6)



### AUROC
data_plot = data_formatted %>% 
  filter(dataset %in% barcoded) %>%
  filter(condition %in% c("classifier", "doublet only control", "singlet only control", "positive control", "shuffled control", "scrambled control", "otherMethods"))
data_plot_summary = data_plot %>% group_by(condition) %>% summarise(sdAuroc = sd(auroc), avgAuroc = mean(auroc))
data_plot <- left_join(data_plot, data_plot_summary, by = "condition")
data_plot$condition <- fct_relevel(data_plot$condition,  "scrambled control", "singlet only control", "doublet only control", "otherMethods", "shuffled control", "positive control", "classifier")

plot = ggplot(data_plot, aes(x = condition, y = auroc)) +
  #geom_violin(data = data_plot, aes(x = condition, y = auroc, fill = condition)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auroc, color = condition), size=1.5, shape = 16, height = 0.01) +
  geom_pointrange(data = data_plot, aes(x = condition, y =avgAuroc, ymin = avgAuroc-sdAuroc, ymax = avgAuroc+sdAuroc, fill = "condition"), color = "black", size = 1) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.04)) + # making sure the error bar fits
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plot_directory, 'classifierControls_AUROC.svg'), width = 8, height = 6)
#ggsave(plot, file = paste0(plot_directory, 'classifierControls_AUROC.png'), width = 8, height = 6)











