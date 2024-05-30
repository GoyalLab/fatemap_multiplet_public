### plotting ensemble method results together for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20240205
### Last edited by Madeline E Melzer on 20240318

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(PRROC)

#BiocManager::install("scds")

set.seed(23)
palette = c("scDblFinder" = "#b26c98", 
            "DoubletFinder" = "#feb325", 
            "Scrublet" = "#2474a9", 
            "hybrid" = "#03bfc4", 
            "scrambled" = "#007024",
            "FM01" = "#2e39ff",
            "FM02" = "#30ebff",
            "ClonMapper" = "#6fa1ff")

plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/ensembleMethods/"

################################
#                              #
#            Chord             #
#                              #
################################

partial_rerun = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/partial_rerun/allBenchmarking_fromPartialRerunFolder_allReranDatasetsUsable.csv")
final = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/final/allBenchmarking_fromFinalFolder_filteredForUsableDatasets.csv")
nonBarcoded = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/final/nonBarcodedDatasets.csv")


# chord, hybrid, and doubletFinder (chord = hybrid + doubletFinder)
# actual doublet rate = 0.08
# datasets: FM01, FM02, and ClonMapper

subfolders <- list.dirs("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/", recursive = FALSE)

chordPlotData <- map_dfr(subfolders, ~{
  detection_files <- list.files(.x, pattern = "all_detection_rates_0.08.tsv", full.names = TRUE)
  # Iterate over each file and read + process the data
  map_dfr(detection_files, ~{
    read_tsv(.x) %>%
      mutate(dataset = ID,
             condition = basename(dirname(.x)),  # Use dirname() to get the folder name
             auroc = roc,
             auprc = pr) %>%
      select(dataset, condition, auroc, auprc) %>%
      filter(str_detect(dataset, "FM01") | str_detect(dataset, "FM02") | str_detect(dataset, "ClonMapper")) %>%
      mutate(condition = str_replace(condition, "doubletFinder", "DoubletFinder")) %>%
      separate(dataset, into = c("dataset", "sample"), sep = "_", extra = "merge")
  })
})

#this file. has FM01 and FM02 data, but the current ClonMapper data is not correct in it- it is pre-re-run. I have filtered out the clonMapper data and added the new one below. 
write_csv(chordPlotData, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chordPlotData.csv")
chordPlotData_old = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chordPlotData.csv")
chordPlotData_old = chordPlotData_old %>% filter(dataset != "ClonMapper")

partial_rerun = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/partial_rerun/allBenchmarking_fromPartialRerunFolder_allReranDatasetsUsable.csv")
benchmarking_ClonMapper = partial_rerun %>% filter(dataset == "ClonMapper") %>% filter(condition %in% c("hybrid", "doublet_finder")) %>% filter(dbl_act == 0.08)  %>%
  mutate(condition = str_replace(condition, "doublet_finder", "DoubletFinder"))


chord_ClonMapper = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chord_ClonMapperToAdd.tsv", sep = "\t")
chord_ClonMapper = chord_ClonMapper %>%
  mutate(auroc = roc,
         auprc = pr) %>%
  separate(ID, into = c("dataset", "sample"), sep = "_", extra = "merge")
chord_ClonMapper$condition = "chord"


chordPlotData = bind_rows(chordPlotData_old, chord_ClonMapper, benchmarking_ClonMapper)
write_csv(chordPlotData, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chordPlotData_new.csv")
chordPlotData = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chordPlotData_new.csv")


chordPlotData_collapsed <- chordPlotData %>%
  group_by(condition, dataset, sample) %>%
  summarise(
    auroc = mean(auroc, na.rm = TRUE),  # Taking the average of auroc, removing NA values if any
    auprc = mean(auprc, na.rm = TRUE)   # Taking the average of auprc, removing NA values if any
  ) %>%
  ungroup()

chordPlotData_summary <- chordPlotData_collapsed %>%
  group_by(condition) %>%
  summarise(
    sdAuprc = sd(auprc, na.rm = TRUE),
    avgAuprc = mean(auprc, na.rm = TRUE)
  )

data_plot = chordPlotData_collapsed

data_plot <- left_join(data_plot, chordPlotData_summary, by = "condition")
mean_data_plot <- data_plot %>%
  group_by(condition, dataset) %>%
  summarise(mean_auprc = mean(auprc, na.rm = TRUE)) %>%
  ungroup()

data_plot$condition = factor(data_plot$condition)
data_plot$condition <- fct_relevel(data_plot$condition, "DoubletFinder", "hybrid", "chord")

plot = ggplot(data_plot, aes(x = condition, y = auprc)) +
  geom_boxplot(data = data_plot, aes(x=condition, y=auprc, color = "#000000"), color = "#000000" , alpha = 0,size=0.5, shape = 16, position = position_dodge(width = 0.3)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auprc, color = dataset),size=1.5, shape = 16, position = position_dodge(width = 0.3)) +
  #geom_line(data = mean_data_plot, aes(x = condition, y = mean_auprc, group = dataset, color = dataset), size = 1, position = position_dodge(width = 0.3)) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot

ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_chord_AUPRC.svg'), width = 3, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_chord_AUPRC.png'), width = 3, height = 4)
  
  
### auroc

chordPlotData_summary <- chordPlotData_collapsed %>%
  group_by(condition) %>%
  summarise(
    sdAuroc = sd(auroc, na.rm = TRUE),  # Change auprc to auroc
    avgAuroc = mean(auroc, na.rm = TRUE)  # Change auprc to auroc
  )

data_plot = chordPlotData_collapsed

data_plot <- left_join(data_plot, chordPlotData_summary, by = "condition")

data_plot$condition = factor(data_plot$condition)
data_plot$condition <- fct_relevel(data_plot$condition, "DoubletFinder", "hybrid", "chord")

plot = ggplot(data_plot, aes(x = condition, y = auroc)) +  # Change auprc to auroc
  geom_boxplot(data = data_plot, aes(x = condition, y = auroc, color = "#000000"), color = "#000000", alpha = 0, size = 0.5, position = position_dodge(width = 0.3)) +  # Change auprc to auroc
  geom_jitter(data = data_plot, aes(x = condition, y = auroc, color = dataset), size = 1.5, shape = 16, position = position_dodge(width = 0.3)) +  # Change auprc to auroc
  scale_color_manual(values = palette) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_chord_AUROC.svg'), width = 3, height = 4)  # Change AUPRC to AUROC in the file name
ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_chord_AUROC.png'), width = 3, height = 4)  # Change AUPRC to AUROC in the file name




################################
#                              #
#            Hybrid            #
#                              #
################################

# plot hybrid, cxds, bcds

### first, extract roc and pr for each method
hybridObjectDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/hybrid/outputObjects/" #already all 0.08 actual doublet rate
rds_files <- list.files(hybridObjectDirectory, pattern = "\\.rds$", full.names = TRUE)
results_df <- data.frame(dataset = character(), sample = character(), method = character(), 
                         roc = numeric(), pr = numeric(), dbl_exp = numeric(), dbl_act = numeric(), stringsAsFactors = FALSE)

for (file_path in rds_files) {
  sce <- readRDS(file_path)
  score <- colData(sce)
  file_name <- basename(file_path)
  parts <- unlist(str_split(file_name, "___|__|\\."))  # Split by '___', '__', and '.'
  dataset <- parts[1]
  sample <- parts[2]
  file_name <- str_remove(file_name, "\\.rds$")
  dbl_exp <- as.numeric(str_extract(file_name, "(?<=exp_)[\\d\\.]+"))
  dbl_act <- as.numeric(str_extract(file_name, "(?<=act_)[\\d\\.]+"))
  
  methods <- c("hybrid", "cxds", "bcds")
  for (method in methods) {
    fg <- score[score$label == "doublet", paste0(method, "_score")]
    bg <- score[score$label == "singlet", paste0(method, "_score")]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
    cur_pr <- pr$auc.integral
    roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
    cur_roc <- roc$auc
    
    # Append results to the dataframe
    results_df <- rbind(results_df, data.frame(dataset = dataset, sample = sample, method = method, 
                                               roc = cur_roc, pr = cur_pr, dbl_exp = dbl_exp, dbl_act = dbl_act))
  }
}

colnames(results_df)[colnames(results_df) == "method"] <- "condition"
colnames(results_df)[colnames(results_df) == "roc"] <- "auroc"
colnames(results_df)[colnames(results_df) == "pr"] <- "auprc"

hybridPlotData = results_df
write_csv(hybridPlotData, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/hybridPlotData_new.csv")
hybridPlotData = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/hybridPlotData_new.csv")
hybridPlotData = hybridPlotData %>% select(-dbl_exp, -dbl_act)
#filteredChordPlotData <- subset(chordPlotData, condition == "DoubletFinder") %>% select(dataset, sample, condition, auprc, auroc)

#hybridPlotData <- rbind(hybridPlotData, filteredChordPlotData)


### auprc

hybridPlotData_collapsed <- hybridPlotData %>%
  group_by(condition, dataset, sample) %>%
  summarise(
    auroc = mean(auroc, na.rm = TRUE),  # Taking the average of auroc, removing NA values if any
    auprc = mean(auprc, na.rm = TRUE)   # Taking the average of auprc, removing NA values if any
  ) %>%
  ungroup()

hybridPlotData_summary <- hybridPlotData_collapsed %>%
  group_by(condition) %>%
  summarise(
    sdAuprc = sd(auprc, na.rm = TRUE),
    avgAuprc = mean(auprc, na.rm = TRUE)
  )

data_plot = hybridPlotData_collapsed
data_plot$condition <- as.factor(data_plot$condition)
data_plot <- left_join(data_plot, hybridPlotData_summary, by = "condition")
mean_data_plot <- data_plot %>%
  group_by(condition, dataset) %>%
  summarise(mean_auprc = mean(auprc, na.rm = TRUE)) %>%
  ungroup()
data_plot$condition <- fct_relevel(data_plot$condition, "bcds", "cxds", "hybrid")
plot = ggplot(data_plot, aes(x = condition, y = auprc)) +
  geom_boxplot(data = data_plot, aes(x=condition, y=auprc, color = "#000000"), color = "#000000" , alpha = 0,size=0.5, shape = 16, position = position_dodge(width = 0.3)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auprc, color = dataset),size=1.5, shape = 16, position = position_dodge(width = 0.3)) +
  #geom_line(data = mean_data_plot, aes(x = condition, y = mean_auprc, group = dataset, color = dataset), size = 1, position = position_dodge(width = 0.3)) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot

ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_hybrid_AUPRC.svg'), width = 3, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_hybrid_AUPRC.png'), width = 3, height = 4)


### auroc

hybridPlotData_summary <- hybridPlotData_collapsed %>%
  group_by(condition) %>%
  summarise(
    sdAuroc = sd(auroc, na.rm = TRUE),  # Change auprc to auroc
    avgAuroc = mean(auroc, na.rm = TRUE)  # Change auprc to auroc
  )

data_plot = hybridPlotData_collapsed
data_plot <- left_join(data_plot, hybridPlotData_summary, by = "condition")

mean_data_plot <- data_plot %>%
  group_by(condition, dataset) %>%
  summarise(mean_auroc = mean(auroc, na.rm = TRUE)) %>%  # Change auprc to auroc
  ungroup()
data_plot$condition <- fct_relevel(data_plot$condition, "bcds", "cxds", "hybrid")
plot = ggplot(data_plot, aes(x = condition, y = auroc)) +  # Change auprc to auroc
  geom_boxplot(data = data_plot, aes(x = condition, y = auroc, color = "#000000"), color = "#000000", alpha = 0, size = 0.5, position = position_dodge(width = 0.3)) +  # Change auprc to auroc
  geom_jitter(data = data_plot, aes(x = condition, y = auroc, color = dataset), size = 1.5, shape = 16, position = position_dodge(width = 0.3)) +  # Change auprc to auroc
  scale_color_manual(values = palette) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_hybrid_AUROC.svg'), width = 3, height = 4)  # Change AUPRC to AUROC in the file name
ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_hybrid_AUROC.png'), width = 3, height = 4)  # Change AUPRC to AUROC in the file name


################################
#                                                          #
# Hybrid and Chord on the same axis (20240318)             #
#                                                          #
################################


chordData = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chordPlotData_new.csv")
hybridData = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/hybridPlotData_new.csv")

chordData = chordData %>% select(dataset, sample, condition, auprc, auroc) %>% filter(condition != "hybrid")
hybridData = hybridData %>% select(dataset, sample, condition, auprc, auroc)

ensembleData = bind_rows(chordData, hybridData)
#write.csv(ensembleData, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chordHybridEnsemblePlotted.csv")


ensembleData = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/ensemblePlotData/chordHybridEnsemblePlotted.csv")


ensembleData_collapsed <- ensembleData %>%
  group_by(condition, dataset, sample) %>%
  summarise(
    auroc = mean(auroc, na.rm = TRUE),  # Taking the average of auroc, removing NA values if any
    auprc = mean(auprc, na.rm = TRUE)   # Taking the average of auprc, removing NA values if any
  ) %>%
  ungroup()


### auprc

ensembleData_summary <- ensembleData_collapsed %>%
  group_by(condition) %>%
  summarise(
    sdAuprc = sd(auprc, na.rm = TRUE),
    avgAuprc = mean(auprc, na.rm = TRUE)
  )

data_plot = ensembleData_collapsed

data_plot <- left_join(data_plot, ensembleData_summary, by = "condition")
mean_data_plot <- data_plot %>%
  group_by(condition, dataset) %>%
  summarise(mean_auprc = mean(auprc, na.rm = TRUE)) %>%
  ungroup()

data_plot$condition = factor(data_plot$condition)
data_plot$condition <- fct_relevel(data_plot$condition, "bcds", "cxds", "hybrid", "DoubletFinder", "chord")

plot = ggplot(data_plot, aes(x = condition, y = auprc)) +
  geom_boxplot(data = data_plot, aes(x=condition, y=auprc, color = "#000000"), color = "#000000" , alpha = 0,size=0.5, shape = 16, position = position_dodge(width = 0.3)) +
  geom_jitter(data = data_plot, aes(x=condition, y=auprc, color = dataset),size=1.5, shape = 16, position = position_dodge(width = 0.3)) +
  #geom_line(data = mean_data_plot, aes(x = condition, y = mean_auprc, group = dataset, color = dataset), size = 1, position = position_dodge(width = 0.3)) +
  scale_color_manual(values = palette) + 
  #scale_fill_manual(values = palette) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        legend.box = "horizontal")
plot

ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_AUPRC.svg'), width = 4.5, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_AUPRC.png'), width = 4.5, height = 3)


### auroc

ensembleData_summary <- ensembleData_collapsed %>%
  group_by(condition) %>%
  summarise(
    sdAuroc = sd(auroc, na.rm = TRUE),  # Change auprc to auroc
    avgAuroc = mean(auroc, na.rm = TRUE)  # Change auprc to auroc
  )

data_plot = ensembleData_collapsed

data_plot <- left_join(data_plot, ensembleData_summary, by = "condition")

data_plot$condition = factor(data_plot$condition)
data_plot$condition <- fct_relevel(data_plot$condition, "bcds", "cxds", "hybrid", "DoubletFinder", "chord")

plot = ggplot(data_plot, aes(x = condition, y = auroc)) +  # Change auprc to auroc
  geom_boxplot(data = data_plot, aes(x = condition, y = auroc, color = "#000000"), color = "#000000", alpha = 0, size = 0.5, position = position_dodge(width = 0.3)) +  # Change auprc to auroc
  geom_jitter(data = data_plot, aes(x = condition, y = auroc, color = dataset), size = 1.5, shape = 16, position = position_dodge(width = 0.3)) +  # Change auprc to auroc
  scale_color_manual(values = palette) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1.01)) +
  labs(y = "", x = "") +
  theme(legend.position = "none")
plot

ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_AUROC.svg'), width = 4.5, height = 3)  # Change AUPRC to AUROC in the file name
ggsave(plot, file = paste0(plotDirectory, 'ensembleMethods_AUROC.png'), width = 4.5, height = 3)  # Change AUPRC to AUROC in the file name



#### statistics 20240321

# using ensembleData_collapsed because that is what is plotted.


#auprc
ensembleStats = ensembleData_collapsed %>% select(condition, auprc, auroc) 

auprc_ensembleStats <- ensembleStats %>%
  select(condition, auprc) %>%
  mutate(row_id = row_number()) %>%
  pivot_wider(names_from = condition, values_from = auprc) %>%
  select(-row_id)

wilcox.test(auprc_ensembleStats$bcds, auprc_ensembleStats$hybrid, paired = FALSE) # W = 6, p-value = 2.219e-05
wilcox.test(auprc_ensembleStats$cxds, auprc_ensembleStats$hybrid, paired = FALSE) # W = 116, p-value = 0.01004
wilcox.test(auprc_ensembleStats$hybrid, auprc_ensembleStats$chord, paired = FALSE) # W = 67, p-value = 0.7987
wilcox.test(auprc_ensembleStats$DoubletFinder, auprc_ensembleStats$chord, paired = FALSE) # W = 74, p-value = 0.9323
wilcox.test(auprc_ensembleStats$bcds, auprc_ensembleStats$chord, paired = FALSE) #W = 1, p-value = 1.479e-06
wilcox.test(auprc_ensembleStats$cxds, auprc_ensembleStats$chord, paired = FALSE) # W = 115, p-value = 0.01209


auroc_ensembleStats <- ensembleStats %>%
  select(condition, auroc) %>%
  mutate(row_id = row_number()) %>%
  pivot_wider(names_from = condition, values_from = auroc) %>%
  select(-row_id)

wilcox.test(auroc_ensembleStats$bcds, auroc_ensembleStats$hybrid, paired = FALSE) # W = 6, p-value = 2.219e-05
wilcox.test(auroc_ensembleStats$cxds, auroc_ensembleStats$hybrid, paired = FALSE) # W = 51, p-value = 0.2415
wilcox.test(auroc_ensembleStats$hybrid, auroc_ensembleStats$chord, paired = FALSE) # W = 54, p-value = 0.3186
wilcox.test(auroc_ensembleStats$DoubletFinder, auroc_ensembleStats$chord, paired = FALSE) # W = 39, p-value = 0.05966
wilcox.test(auroc_ensembleStats$bcds, auroc_ensembleStats$chord, paired = FALSE) # W = 5, p-value = 1.405e-05
wilcox.test(auroc_ensembleStats$cxds, auroc_ensembleStats$chord, paired = FALSE) # W = 38, p-value = 0.05186




