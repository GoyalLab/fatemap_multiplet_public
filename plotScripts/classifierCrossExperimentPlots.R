### plotting cross-experimental classifier results
### Created by Madeline E Melzer on 20240513
### Last edited by Madeline E Melzer on 20240515

library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(ggbreak)
library(scales)

set.seed(23)

plot_directory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/classifier/crossExperiment/"
dataset_dir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/"


palette = c("classifier" = "#f59fff", 
            "scDblFinder" = "#b26c98", 
            "DoubletFinder" = "#feb325", 
            "Scrublet" = "#2474a9", 
            "hybrid" = "#03bfc4", 
            "scrambled" = "#007024",
            "otherMethods" = "#808285",
            "sample A classifier" = "#006838",
            "sample B classifier" = "#8dc63f")

# load doublet detection method results for comparison, formatting properly
ddmResults = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/finalBenchmarking/reDownloaded/barcodedAllForPlot.csv") %>% 
  filter(dbl_act == 0.1) %>%
  dplyr::rename(technology = dataset) %>%
  select(-X) %>%
  mutate(classifier = FALSE)

# load individual classifier results
dataset = "ClonMapper"
clonMapper = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/ClonMapper/summary_s1s2_bothClassifiers.csv")
clonMapper = clonMapper %>% mutate(technology = "ClonMapper") %>% mutate(classifier = TRUE)
ddmClonMapper = ddmResults %>% filter(technology == "ClonMapper") %>% filter(sample %in% c("FM1", "FM7")) %>%
  mutate(sample = case_when(
    sample == "FM1" ~ "sample1",
    sample == "FM7" ~ "sample2",
    TRUE ~ sample  # This line keeps other sample values unchanged if needed
  ))
clonMapper = bind_rows(clonMapper, ddmClonMapper)

dataset = "TREX"
trex = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/TREX/summary_s1s2_bothClassifiers.csv")
trex_s1 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/TREX/summary_s2_s1Classifier_test4.csv")
trex_s2 = filter(trex, condition == "sample2")
trex = bind_rows(trex_s1, trex_s2)
trex = trex %>% mutate(technology = "TREX") %>% mutate(classifier = TRUE)
ddmTREX = ddmResults %>% filter(technology == "TREX") %>% filter(sample %in% c("brain2", "brain3")) %>%
  mutate(sample = case_when(
    sample == "brain2" ~ "sample1",
    sample == "brain3" ~ "sample2",
    TRUE ~ sample  # This line keeps other sample values unchanged if needed
  ))
trex = bind_rows(trex, ddmTREX)

dataset = "TREX_minusCluster"
trex_minusCluster = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/TREX_minusCluster/summary_s1s2_bothClassifiers.csv")
trex_minusCluster = trex_minusCluster %>% mutate(technology = "TREX_minusCluster") %>% mutate(classifier = TRUE)
ddmTREX = ddmResults %>% filter(technology == "TREX") %>% filter(sample %in% c("brain2", "brain3")) %>%
  mutate(sample = case_when(
    sample == "brain2" ~ "sample1",
    sample == "brain3" ~ "sample2",
    TRUE ~ sample  # This line keeps other sample values unchanged if needed
  ))
trex_minusCluster = bind_rows(trex_minusCluster, ddmTREX)

dataset = "SPLINTR"
splintr_veh = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/SPLINTR/summary_s1s2_bothClassifiers.csv")
splintr_veh = splintr_veh %>% mutate(technology = "SPLINTR_veh") %>% mutate(classifier = TRUE)
ddmSPLINTR_veh = ddmResults %>% filter(technology == "SPLINTR") %>% filter(sample %in% c("chemoVehicle_1", "chemoVehicle_2")) %>%
  mutate(technology = "SPLINTR_veh",
    sample = case_when(
    sample == "chemoVehicle_1" ~ "sample1",
    sample == "chemoVehicle_2" ~ "sample2",
    TRUE ~ sample  # This line keeps other sample values unchanged if needed
  ))
splintr_veh = bind_rows(splintr_veh, ddmSPLINTR_veh)

splintr_chemo = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/SPLINTR/summary_s3s4_bothClassifiers.csv")
splintr_chemo = splintr_chemo %>% mutate(classifier = TRUE) %>%
  mutate(technology = "SPLINTR_chemo",
    dataset = case_when(
      dataset == "sample3" ~ "sample1", #sample 3 in overall SPLINTR, but doing this for coloring
      dataset == "sample4" ~ "sample2", #sample 4 in overall SPLINTR, but doing this for coloring
    TRUE ~ dataset),
    condition = case_when(
      condition == "sample3" ~ "sample1", #sample 3 in overall SPLINTR, but doing this for coloring
      condition == "sample4" ~ "sample2", #sample 4 in overall SPLINTR, but doing this for coloring
    TRUE ~ condition  # This line keeps other sample values unchanged if needed
  ))
ddmSPLINTR_chemo = ddmResults %>% filter(technology == "SPLINTR") %>% filter(sample %in% c("chemoDay2_1", "chemoDay2_2")) %>%
  mutate(technology = "SPLINTR_chemo",
    sample = case_when(
    sample == "chemoDay2_1" ~ "sample1", #sample 3 in overall SPLINTR, but doing this for coloring
    sample == "chemoDay2_2" ~ "sample2", #sample 4 in overall SPLINTR, but doing this for coloring
    TRUE ~ sample  # This line keeps other sample values unchanged if needed
  ))
splintr_chemo = bind_rows(splintr_chemo, ddmSPLINTR_chemo)

dataset = "Goyal et al. 1"
fateMap_data = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/singletCode20240326/plotData/classifier/alls1s2ForPlot_formatted.csv")
fateMap_data = fateMap_data %>% mutate(row_number = row_number(),  technology = "Goyal et al. 1", classifier = if_else(condition %in% c("sample1", "sample2"), TRUE, FALSE)) #adding row number so i can know what classifier it goes to. 
set.seed(23) #need to set seed IMMEDIATELY BEFORE sample_n() otherwise it will still be random. 
downsampled_fateMap <- fateMap_data %>% filter(condition %in% c("sample1", "sample2")) %>% group_by(condition) %>% sample_n(10) %>% ungroup()
#getting classifier numbers to run them on scrambled sample controls. 
downsampled_fateMap_numbered = downsampled_fateMap %>% mutate(classifierNumber = row_number - 1, classifierNumber = if_else(condition == "sample1", classifierNumber - 100, classifierNumber))
classifierNumbers_sample1 <- downsampled_fateMap_numbered %>% filter(condition == "sample1") %>% pull(classifierNumber)
classifierNumbers_sample2 <- downsampled_fateMap_numbered %>% filter(condition == "sample2") %>% pull(classifierNumber)
print(classifierNumbers_sample1) #[1] 28 27 71 42 44 33 47 16 20 92
print(classifierNumbers_sample2) # [1] 39 35 69 97 30 85 76 41  9 88
ddmFateMap = fateMap_data %>% filter(condition %in% c("DoubletFinder", "scDblFinder", "hybrid", "Scrublet"))
fateMap = bind_rows(downsampled_fateMap, ddmFateMap)


######## combining all classifiers for cell type plots

allSeparatePairs = bind_rows(clonMapper, trex, splintr_veh, splintr_chemo, fateMap)

allSeparatePairs_formatted <- allSeparatePairs %>%
  mutate(condition = str_replace(condition, ".*neg_control_doublets.*", "doublet control"),
         condition = str_replace(condition, ".*neg_control_singlets.*", "singlet control"),
         condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
         condition = str_replace(condition, ".*shuffled.*", "shuffled"),
         #condition = if_else(condition == as.character(dataset), "classifier", condition))
         condition = str_replace(condition, "sample1", "sample A classifier"),
         condition = str_replace(condition, "sample2", "sample B classifier"),
         dataset = if_else(is.na(dataset), sample, dataset)) %>%
  mutate(technology_dataset = paste(technology, dataset, sep = "_")) # making dataset-technology pair


# plot of each dataset and classifier, no classifier
ggplot(allSeparatePairs_formatted, aes(x = technology, y = auprc, color = condition, group = interaction(technology, classifier))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  # Draw boxplots
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 2, shape = 16) +  # Add jittered points
  scale_color_manual(values = palette) +
  labs( y = "AUPRC",
        color = "method") +
  theme_classic()

# all datasets combined classifier, no classifier
classifier_labels <- c(`TRUE` = "classifier", `FALSE` = "other methods")

ggplot(allSeparatePairs_formatted, aes(x = classifier, y = auprc, color = condition, group = as.factor(classifier))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  # Draw boxplots
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 2, shape = 16) +  # Add jittered points
  scale_color_manual(values = palette) +
  scale_x_discrete(labels = classifier_labels) +
  labs( y = "AUPRC",
        x = "",
        color = "method") +
  theme_classic()

# statistics for all datasets combined classifier vs. other methods
auprc_classifier <- allSeparatePairs_formatted %>%
  filter(classifier == TRUE) %>%
  pull(auprc)
auprc_otherMethods <- allSeparatePairs_formatted %>%
  filter(classifier == FALSE) %>%
  pull(auprc)

# Perform the Wilcoxon rank-sum test
wilcox.test(auprc_classifier, auprc_otherMethods, paired = FALSE) #W = 11070, p-value = 9.81e-09, 20240514
t.test(auprc_classifier, auprc_otherMethods, paired = FALSE) #t = 6.8287, df = 206.1, p-value = 9.426e-11, 20240514


### same cell type, different conditions with SPLINTR vehicle vs chemo day 2 (did not use)

splintr_crossCondition = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/SPLINTR/summary_s1s3_bothClassifiers.csv")
splintr_crossCondition = splintr_crossCondition %>% mutate(classifier = TRUE) %>%
  mutate(technology = "SPLINTR_crossCondition",
         dataset = case_when(
           dataset == "sample1_s1s2s3s4" ~ "sample1", #sample 3 in overall SPLINTR, but doing this for coloring
           dataset == "sample3_s1s2s3s4" ~ "sample2", #sample 4 in overall SPLINTR, but doing this for coloring
           TRUE ~ dataset),
         condition = case_when(
           condition == "sample1_s1s2s3s4" ~ "sample1", #sample 3 in overall SPLINTR, but doing this for coloring
           condition == "sample3_s1s2s3s4" ~ "sample2", #sample 4 in overall SPLINTR, but doing this for coloring
           TRUE ~ condition  # This line keeps other sample values unchanged if needed
         ))
ddmSPLINTR_crossCondition = ddmResults %>% filter(technology == "SPLINTR") %>% filter(sample %in% c("chemoVehicle_1", "chemoDay2_1")) %>%
  mutate(technology = "SPLINTR_crossCondition",
         sample = case_when(
           sample == "chemoVehicle_1" ~ "sample1", #sample 3 in overall SPLINTR, but doing this for coloring
           sample == "chemoDay2_1" ~ "sample2", #sample 4 in overall SPLINTR, but doing this for coloring
           TRUE ~ sample  # This line keeps other sample values unchanged if needed
         ))

splintr_crossCondition = bind_rows(splintr_crossCondition, ddmSPLINTR_crossCondition)

crossCondition_formatted <- splintr_crossCondition %>%
  mutate(condition = str_replace(condition, "sample1", "sample A classifier"),
         condition = str_replace(condition, "sample2", "sample B classifier"),
         dataset = if_else(is.na(dataset), sample, dataset)) # making dataset-technology pair

# plot of each dataset and classifier, no classifier
ggplot(crossCondition_formatted, aes(x = technology, y = auprc, color = condition, group = interaction(technology, classifier))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  # Draw boxplots
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 2, shape = 16) +  # Add jittered points
  scale_color_manual(values = palette) +
  labs( y = "AUPRC",
        x = "",
        color = "method") +
  theme_classic()




### different cell type, conditions, technology in humans with FateMap and ClonMapper samples

human_allA = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/summary_CM2FM2_humanAClassifier.csv")
human_allB = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/summary_CM1FM1_humanBClassifier.csv")
human_all = bind_rows(human_allA, human_allB)
human_all = human_all %>% mutate(classifier = TRUE) %>%
  mutate(technology = "human_all")
ddmHumanAll = ddmResults %>% filter(technology %in% c("FM01", "ClonMapper")) %>% filter(sample %in% c("sample1", "FM1", "sample2", "FM7")) %>%
  mutate(technology = "human_all")
human_all = bind_rows(human_all, ddmHumanAll)

human_all_formatted <- human_all %>%
  mutate(dataset = if_else(is.na(dataset), sample, dataset)) # making dataset-technology pair

# plot of each dataset and classifier, no classifier
ggplot(human_all_formatted, aes(x = technology, y = auprc, color = condition, group = interaction(technology, classifier))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  # Draw boxplots
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 2, shape = 16) +  # Add jittered points
  scale_color_manual(values = palette) +
  labs( y = "AUPRC",
        x = "",
        color = "method") +
  theme_classic()


### different cell type, conditions, technology in mice with SPLINTR and TREX samples

mouse_all1 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/summary_SP2SP4TX2_mouseAClassifier.csv")
mouse_all2 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/summary_SP1SP3TX1_mouseBClassifier.csv")
mouse_all = bind_rows(mouse_all1, mouse_all2)
mouse_all = mouse_all %>% mutate(classifier = TRUE) %>%
  mutate(technology = "mouse_all")
ddmMouseAll = ddmResults %>% filter(technology %in% c("SPLINTR", "TREX")) %>% filter(sample %in% c("chemoVehicle_1", "chemoVehicle_2", "chemoDay2_1", "chemoDay2_2", "brain2", "brain3")) %>%
  mutate(technology = "mouse_all")
mouse_all = bind_rows(mouse_all, ddmMouseAll)

mouse_all_formatted <- mouse_all %>%
  mutate(dataset = if_else(is.na(dataset), sample, dataset)) # making dataset-technology pair

# plot of each dataset and classifier, no classifier
ggplot(mouse_all_formatted, aes(x = technology, y = auprc, color = condition, group = interaction(technology, classifier))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  # Draw boxplots
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 2, shape = 16) +  # Add jittered points
  scale_color_manual(values = palette) +
  labs( y = "AUPRC",
        x = "",
        color = "method") +
  theme_classic()


### plotting all different classifier tests on one axis

allClassifiers = bind_rows(allSeparatePairs_formatted, human_all_formatted, mouse_all_formatted) %>%
  mutate(technology = factor(technology, levels = c("Goyal et al. 1", "ClonMapper", "TREX", 
                                                    "SPLINTR_veh", "SPLINTR_chemo", 
                                                    "SPLINTR_crossCondition", "human_all", "mouse_all")))
#write.csv(allClassifiers, paste0(dataset_dir, "allCrossExperimentClassifiersAndDDMResults.csv")) #20240517 MEM

#allClassifiers_summary = allClassifiers %>% group_by(technology, classifier) %>% summarise(sdAuprc = sd(auprc), avgAuprc = mean(auprc)) #only need these if using geom_pointrange
#allClassifiers_plot <- left_join(allClassifiers, allClassifiers_summary, by = c("technology", "classifier"))

ggplot(allClassifiers, aes(x = technology, y = auprc, color = condition, group = interaction(technology, classifier))) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 2, shape = 16) +  
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUPRC",
        x = "",
        color = "method") +
  theme_classic()



#### plotting cell-type differences figure

cellTypeSpecificClassifiers = allClassifiers %>% filter(technology %in% c("Goyal et al. 1", "ClonMapper", "TREX", "SPLINTR_veh", "SPLINTR_chemo")) %>%
  mutate(technology = case_when(
    technology == "SPLINTR_veh" ~ "SPLINTR",
    technology == "SPLINTR_chemo" ~ "SPLINTR",
    TRUE ~ technology))  %>%
  mutate(technology = factor(technology, levels = c("Goyal et al. 1", "ClonMapper", "TREX", "SPLINTR")))

# AUPRC
plot = ggplot(cellTypeSpecificClassifiers, aes(x = technology, y = auprc, color = condition, group = interaction(technology, classifier))) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUPRC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'cellTypeSpecificClassifiers_auprc.svg'), width = 5, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'cellTypeSpecificClassifiers_auprc.png'), width = 5, height = 4) #20240517 MEM

# AUROC
plot = ggplot(cellTypeSpecificClassifiers, aes(x = technology, y = auroc, color = condition, group = interaction(technology, classifier))) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUROC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'cellTypeSpecificClassifiers_auroc.svg'), width = 5, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'cellTypeSpecificClassifiers_auroc.png'), width = 5, height = 4) #20240517 MEM



#### plotting all human and all mouse classifier figure

allHumanMouseClassifiers = allClassifiers %>% filter(technology %in% c("human_all", "mouse_all"))

# AUPRC
plot = ggplot(allHumanMouseClassifiers, aes(x = technology, y = auprc, color = condition, group = interaction(technology, classifier))) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUPRC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auprc.svg'), width = 3, height = 4) #20240515 MEM
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auprc.png'), width = 3, height = 4) #20240515 MEM

# AUROC
plot = ggplot(allHumanMouseClassifiers, aes(x = technology, y = auroc, color = condition, group = interaction(technology, classifier))) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUROC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auroc.svg'), width = 3, height = 4) #20240515 MEM
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auroc.png'), width = 3, height = 4) #20240515 MEM



############## statistics (Wilcox and t-test but only using Wilcox)


### AUPRC

# Goyal et al. 1 (note: different than 6H because downsampled to n = 10 to match other dataset n)
auprc_classifier_FateMap <- allClassifiers %>% filter(technology == "Goyal et al. 1") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_FateMap <- allClassifiers %>% filter(technology == "Goyal et al. 1") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_FateMap, auprc_otherMethods_FateMap, paired = FALSE) #W = 959, p-value = 1.189e-10 #20240517 MEM
t.test(auprc_classifier_FateMap, auprc_otherMethods_FateMap, paired = FALSE) #t = 8.8757, df = 21.728, p-value = 1.122e-08, #20240517 MEM

# ClonMapper
auprc_classifier_ClonMapper <- allClassifiers %>% filter(technology == "ClonMapper") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_ClonMapper <- allClassifiers %>% filter(technology == "ClonMapper") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_ClonMapper, auprc_otherMethods_ClonMapper, paired = FALSE) #W = 960, p-value = 1.081e-10 #20240515 MEM
t.test(auprc_classifier_ClonMapper, auprc_otherMethods_ClonMapper, paired = FALSE) #t = 12.816, df = 42.703, p-value = 3.138e-16, #20240515 MEM

# TREX
auprc_classifier_TREX <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_TREX <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_TREX, auprc_otherMethods_TREX, paired = FALSE) #W = 304, p-value = 0.01808 #20240517 MEM
t.test(auprc_classifier_TREX, auprc_otherMethods_TREX, paired = FALSE) #t = -2.8125, df = 23.964, p-value = 0.009652, #20240517 MEM

auprc_classifier_TREX_s1 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == TRUE) %>% filter(condition == "sample A classifier") %>% pull(auprc)
auprc_otherMethods_TREX_s1 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == FALSE) %>% filter(sample == "sample1") %>% pull(auprc)

wilcox.test(auprc_classifier_TREX_s1, auprc_otherMethods_TREX_s1, paired = FALSE) #W = 159, p-value = 0.1447 #20240517 MEM
t.test(auprc_classifier_TREX_s1, auprc_otherMethods_TREX_s1, paired = FALSE) #t = 2.3071, df = 31.795, p-value = 0.02772, #20240517 MEM

auprc_classifier_TREX_s2 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == TRUE) %>% filter(condition == "sample B classifier") %>% pull(auprc)
auprc_otherMethods_TREX_s2 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == FALSE) %>% filter(sample == "sample2") %>% pull(auprc)

wilcox.test(auprc_classifier_TREX_s2, auprc_otherMethods_TREX_s2, paired = FALSE) #W = 0, p-value = 5.825e-06 #20240517 MEM
t.test(auprc_classifier_TREX_s2, auprc_otherMethods_TREX_s2, paired = FALSE) #t = -11.362, df = 17.806, p-value = 1.368e-09, #20240517 MEM

# SPLINTR_veh
auprc_classifier_SPLINTR_veh <- allClassifiers %>% filter(technology == "SPLINTR_veh") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_SPLINTR_veh <- allClassifiers %>% filter(technology == "SPLINTR_veh") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_SPLINTR_veh, auprc_otherMethods_SPLINTR_veh, paired = FALSE) #W = 763, p-value = 0.0001431 #20240515 MEM
t.test(auprc_classifier_SPLINTR_veh, auprc_otherMethods_SPLINTR_veh, paired = FALSE) #t = 5.2671, df = 61.141, p-value = 1.899e-06, #20240515 MEM

# SPLINTR_chemo
auprc_classifier_SPLINTR_chemo <- allClassifiers %>% filter(technology == "SPLINTR_chemo") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_SPLINTR_chemo <- allClassifiers %>% filter(technology == "SPLINTR_chemo") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_SPLINTR_chemo, auprc_otherMethods_SPLINTR_chemo, paired = FALSE) #W = 861, p-value = 3.027e-07 #20240515 MEM
t.test(auprc_classifier_SPLINTR_chemo, auprc_otherMethods_SPLINTR_chemo, paired = FALSE) #t = 6.6164, df = 52.551, p-value = 1.936e-08, #20240515 MEM

# SPLINTR (both veh and chemo)
auprc_classifier_SPLINTR <- allClassifiers %>% filter(technology %in% c("SPLINTR_chemo", "SPLINTR_veh")) %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_SPLINTR <- allClassifiers %>% filter(technology == c("SPLINTR_chemo", "SPLINTR_veh")) %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_SPLINTR, auprc_otherMethods_SPLINTR, paired = FALSE) #W = 1466, p-value = 2.273e-05 #20240515 MEM
t.test(auprc_classifier_SPLINTR, auprc_otherMethods_SPLINTR, paired = FALSE) #t = 5.5742, df = 79.696, p-value = 3.279e-07, #20240515 MEM

# SPLINTR_crossCondition
auprc_classifier_SPLINTR_crossCondition <- allClassifiers %>% filter(technology == "SPLINTR_crossCondition") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_SPLINTR_crossCondition <- allClassifiers %>% filter(technology == "SPLINTR_crossCondition") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_SPLINTR_crossCondition, auprc_otherMethods_SPLINTR_crossCondition, paired = FALSE) #W = 820, p-value = 4.879e-06 #20240515 MEM
t.test(auprc_classifier_SPLINTR_crossCondition, auprc_otherMethods_SPLINTR_crossCondition, paired = FALSE) #t = 5.3331, df = 64.562, p-value = 1.323e-06, #20240515 MEM

# human_all
auprc_classifier_human_all <- allClassifiers %>% filter(technology == "human_all") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_human_all <- allClassifiers %>% filter(technology == "human_all") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_human_all, auprc_otherMethods_human_all, paired = FALSE) #W = 3792, p-value < 2.2e-16 #20240515 MEM
t.test(auprc_classifier_human_all, auprc_otherMethods_human_all, paired = FALSE) #t = 12.515, df = 66.075, p-value < 2.2e-16, #20240515 MEM

mean(auprc_classifier_human_all) #0.575315 #20240515 MEM
mean(auprc_otherMethods_human_all) #0.3228827 #20240515 MEM

# mouse_all
auprc_classifier_mouse_all <- allClassifiers %>% filter(technology == "mouse_all") %>% filter(classifier == TRUE) %>% pull(auprc)
auprc_otherMethods_mouse_all <- allClassifiers %>% filter(technology == "mouse_all") %>% filter(classifier == FALSE) %>% pull(auprc)

wilcox.test(auprc_classifier_mouse_all, auprc_otherMethods_mouse_all, paired = FALSE) #W = 5425, p-value = 0.00404 #20240515 MEM
t.test(auprc_classifier_mouse_all, auprc_otherMethods_mouse_all, paired = FALSE) #t = 4.1982, df = 163.01, p-value = 4.413e-05, #20240515 MEM

mean(auprc_classifier_mouse_all) #0.5701359 #20240515 MEM
mean(auprc_otherMethods_mouse_all) #0.4686916 #20240515 MEM


# testing to see if classifier efficacy changes when you add more conditions, or when you add everything together.

wilcox.test(auprc_classifier_SPLINTR_veh, auprc_classifier_SPLINTR_crossCondition, paired = FALSE) #W = 232, p-value = 0.3983 #20240515 MEM
wilcox.test(auprc_classifier_SPLINTR_chemo, auprc_classifier_SPLINTR_crossCondition, paired = FALSE) #W = 272, p-value = 0.0524 #20240515 MEM

wilcox.test(auprc_classifier_FateMap, auprc_classifier_human_all, paired = FALSE) #W = 612, p-value = 0.00066 #20240515 MEM
wilcox.test(auprc_classifier_ClonMapper, auprc_classifier_human_all, paired = FALSE) #W = 288, W = 532, p-value = 0.03845 #20240515 MEM


### AUROC

# Goyal et al. 1 (note: different than 6H because downsampled to n = 10 to match other dataset n)
auroc_classifier_FateMap <- allClassifiers %>% filter(technology == "Goyal et al. 1") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_FateMap <- allClassifiers %>% filter(technology == "Goyal et al. 1") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_FateMap, auroc_otherMethods_FateMap, paired = FALSE) #W = 960, p-value = 1.088e-10 #20240517 MEM
t.test(auroc_classifier_FateMap, auroc_otherMethods_FateMap, paired = FALSE) #t = 12.234, df = 35.018, p-value = 3.351e-14, #20240517 MEM

# ClonMapper
auroc_classifier_ClonMapper <- allClassifiers %>% filter(technology == "ClonMapper") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_ClonMapper <- allClassifiers %>% filter(technology == "ClonMapper") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_ClonMapper, auroc_otherMethods_ClonMapper, paired = FALSE) #W = 959, p-value = 1.181e-10 #20240515 MEM
t.test(auroc_classifier_ClonMapper, auroc_otherMethods_ClonMapper, paired = FALSE) #t = 9.5749, df = 61.061, p-value = 9.052e-14, #20240515 MEM

# TREX
auroc_classifier_TREX <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_TREX <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_TREX, auroc_otherMethods_TREX, paired = FALSE) #W = 423, p-value = 0.4466 #20240517 MEM
t.test(auroc_classifier_TREX, auroc_otherMethods_TREX, paired = FALSE) #t = -1.7989, df = 24.194, p-value = 0.08452, #20240517 MEM

auroc_classifier_TREX_s1 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == TRUE) %>% filter(condition == "sample A classifier") %>% pull(auroc)
auroc_otherMethods_TREX_s1 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == FALSE) %>% filter(sample == "sample1") %>% pull(auroc)

wilcox.test(auroc_classifier_TREX_s1, auroc_otherMethods_TREX_s1, paired = FALSE) #W = 205, p-value = 0.001368 #20240517 MEM
t.test(auroc_classifier_TREX_s1, auroc_otherMethods_TREX_s1, paired = FALSE) #t = 3.5897, df = 28.901, p-value = 0.001207, #20240517 MEM

auroc_classifier_TREX_s2 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == TRUE) %>% filter(condition == "sample B classifier") %>% pull(auroc)
auroc_otherMethods_TREX_s2 <- allClassifiers %>% filter(technology == "TREX") %>% filter(classifier == FALSE) %>% filter(sample == "sample2") %>% pull(auroc)

wilcox.test(auroc_classifier_TREX_s2, auroc_otherMethods_TREX_s2, paired = FALSE) #W = 4, p-value = 1.182e-05 #20240517 MEM
t.test(auroc_classifier_TREX_s2, auroc_otherMethods_TREX_s2, paired = FALSE) #t = -8.1517, df = 16.403, p-value = 3.642e-07, #20240517 MEM

# SPLINTR_veh
auroc_classifier_SPLINTR_veh <- allClassifiers %>% filter(technology == "SPLINTR_veh") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_SPLINTR_veh <- allClassifiers %>% filter(technology == "SPLINTR_veh") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_SPLINTR_veh, auroc_otherMethods_SPLINTR_veh, paired = FALSE) #W = 719, p-value = 0.001325 #20240515 MEM
t.test(auroc_classifier_SPLINTR_veh, auroc_otherMethods_SPLINTR_veh, paired = FALSE) #t = 4.6304, df = 51.267, p-value = 2.527e-05, #20240515 MEM

# SPLINTR_chemo
auroc_classifier_SPLINTR_chemo <- allClassifiers %>% filter(technology == "SPLINTR_chemo") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_SPLINTR_chemo <- allClassifiers %>% filter(technology == "SPLINTR_chemo") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_SPLINTR_chemo, auroc_otherMethods_SPLINTR_chemo, paired = FALSE) #W = 947, p-value = 3.4e-10 #20240515 MEM
t.test(auroc_classifier_SPLINTR_chemo, auroc_otherMethods_SPLINTR_chemo, paired = FALSE) #t = 6.9631, df = 51.211, p-value = 6.101e-09, #20240515 MEM

# SPLINTR (both veh and chemo)
auroc_classifier_SPLINTR <- allClassifiers %>% filter(technology %in% c("SPLINTR_chemo", "SPLINTR_veh")) %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_SPLINTR <- allClassifiers %>% filter(technology == c("SPLINTR_chemo", "SPLINTR_veh")) %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_SPLINTR, auroc_otherMethods_SPLINTR, paired = FALSE) #W = 1587, p-value = 1.519e-07 #20240515 MEM
t.test(auroc_classifier_SPLINTR, auroc_otherMethods_SPLINTR, paired = FALSE) #t = 5.2166, df = 50.556, p-value = 3.419e-06, #20240515 MEM

# SPLINTR_crossCondition
auroc_classifier_SPLINTR_crossCondition <- allClassifiers %>% filter(technology == "SPLINTR_crossCondition") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_SPLINTR_crossCondition <- allClassifiers %>% filter(technology == "SPLINTR_crossCondition") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_SPLINTR_crossCondition, auroc_otherMethods_SPLINTR_crossCondition, paired = FALSE) #W = 719, p-value = 0.001326 #20240515 MEM
t.test(auroc_classifier_SPLINTR_crossCondition, auroc_otherMethods_SPLINTR_crossCondition, paired = FALSE) #t = 4.4717, df = 51.878, p-value = 4.251e-05, #20240515 MEM

# human_all
auroc_classifier_human_all <- allClassifiers %>% filter(technology == "human_all") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_human_all <- allClassifiers %>% filter(technology == "human_all") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_human_all, auroc_otherMethods_human_all, paired = FALSE) #W = 3840, p-value < 2.2e-16 #20240515 MEM
t.test(auroc_classifier_human_all, auroc_otherMethods_human_all, paired = FALSE) #t = 13.187, df = 131.78, p-value < 2.2e-16, #20240515 MEM

mean(auroc_classifier_human_all) #0.9039119 #20240515 MEM
mean(auroc_otherMethods_human_all) #0.7592817 #20240515 MEM

# mouse_all
auroc_classifier_mouse_all <- allClassifiers %>% filter(technology == "mouse_all") %>% filter(classifier == TRUE) %>% pull(auroc)
auroc_otherMethods_mouse_all <- allClassifiers %>% filter(technology == "mouse_all") %>% filter(classifier == FALSE) %>% pull(auroc)

wilcox.test(auroc_classifier_mouse_all, auroc_otherMethods_mouse_all, paired = FALSE) #W = 5807, p-value = 0.0001091 #20240515 MEM
t.test(auroc_classifier_mouse_all, auroc_otherMethods_mouse_all, paired = FALSE) #t = 6.0438, df = 175.95, p-value = 8.766e-09, #20240515 MEM

mean(auroc_classifier_mouse_all) #0.8926778 #20240515 MEM
mean(auroc_otherMethods_mouse_all) #0.799303 #20240515 MEM

# testing to see if classifier efficacy changes when you add more conditions, or when you add everything together.

wilcox.test(auroc_classifier_SPLINTR_veh, auroc_classifier_SPLINTR_crossCondition, paired = FALSE) #W = 282, p-value = 0.02633 #20240515 MEM
wilcox.test(auroc_classifier_SPLINTR_chemo, auroc_classifier_SPLINTR_crossCondition, paired = FALSE) #W = 383, p-value = 1.758e-08 #20240515 MEM

wilcox.test(auroc_classifier_FateMap, auroc_classifier_human_all, paired = FALSE) #W = 513, p-value = 0.07761 #20240515 MEM
wilcox.test(auroc_classifier_ClonMapper, auroc_classifier_human_all, paired = FALSE) #W = 558, p-value = 0.01263 #20240515 MEM









##################### classifiers built with (SCRAMBLED CLASSIFIER) or tested on (SCRAMBLED TEST SAMPLE) scrambled samples as a negative control

### for cell type specific classifiers ###

# SCRAMBLED CLASSIFIER

dataset = "ClonMapper"
sample1_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample1.csv"))
sample2_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample2.csv"))
clonMapper_scrambClass = bind_rows(sample1_scrambClass, sample2_scrambClass)

dataset = "TREX"
sample1_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample1.csv"))
sample2_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample2.csv"))
trex_scrambClass = bind_rows(sample1_scrambClass, sample2_scrambClass)

dataset = "SPLINTR"
sample1_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample1.csv"))
sample2_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample2.csv"))
sample3_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample3.csv"))
sample4_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample4.csv"))
splintr_scrambClass = bind_rows(sample1_scrambClass, sample2_scrambClass, sample3_scrambClass, sample4_scrambClass)

dataset = "FateMap" #Goyal et al. 1
s1s2_scrambClass = read.csv(paste0(dataset_dir, dataset, "/sample1_sample2_classifiers.csv"), header = FALSE, stringsAsFactors = FALSE)
s1s2_scrambClass = s1s2_scrambClass %>% select(V1, V2, V3, V4, V5)
colnames(s1s2_scrambClass) = as.character(s1s2_scrambClass[1, ])
s1s2_scrambClass = s1s2_scrambClass[-1, ]
s1s2_scrambClass_trimmed = s1s2_scrambClass[-(1:100), ]
s1s2_scrambClass_trimmed = s1s2_scrambClass_trimmed %>% filter(dataset %in% c("sample1", "sample2")) %>% arrange(dataset)
row.names(s1s2_scrambClass_trimmed) <- NULL
#need to select the correct corresponding classifiers that were randomly chosen above (from fateMap_data)
s1s2_scrambClass_trimmed$classifierNumber <- rep(1:(nrow(s1s2_scrambClass_trimmed)/5), each = 5)
s1s2_scrambClass_trimmed = s1s2_scrambClass_trimmed %>% mutate(classifierNumber = if_else(dataset == "sample2", classifierNumber - 100, classifierNumber))
sample1 <- s1s2_scrambClass_trimmed %>%
  filter(dataset == "sample1" & classifierNumber %in% classifierNumbers_sample1) #28 27 71 42 44 33 47 16 20 92
sample2 <- s1s2_scrambClass_trimmed %>%
  filter(dataset == "sample2" & classifierNumber %in% classifierNumbers_sample2) #39 35 69 97 30 85 76 41  9 88
fateMap_scrambClass = bind_rows(sample1, sample2) %>% mutate(auroc = as.numeric(auroc), auprc = as.numeric(auprc), accuracy = as.numeric(accuracy))



scrambledClassifier_cellType = bind_rows(clonMapper_scrambClass, trex_scrambClass, splintr_scrambClass, fateMap_scrambClass) %>% mutate(condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
                                                                                                                                  condition = if_else(condition == as.character(dataset), "classifier", condition))
scrambledClassifier_cellType = scrambledClassifier_cellType  %>% filter(condition %in% c("scrambled", "classifier"))

#write.csv(scrambledClassifier_cellType, file = paste0(dataset_dir, "allCellTypeClassifiersScrambledClassifierResults.csv")) #20240517 MEM

# AUPRC
plot = ggplot(scrambledClassifier_cellType, aes(x = condition, y = auprc, color = condition)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUPRC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auprc_scrambledClassifier.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auprc_scrambledClassifier.png'), width = 2, height = 4) #20240517 MEM

# AUROC
plot = ggplot(scrambledClassifier_cellType, aes(x = condition, y = auroc, color = condition)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUROC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auroc_scrambledClassifier.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auroc_scrambledClassifier.png'), width = 2, height = 4) #20240517 MEM

### STATISTICS

# AUPRC
auprc_classifier <- scrambledClassifier_cellType %>% filter(condition == "classifier") %>% pull(auprc)
auprc_scrambledClassifier <- scrambledClassifier_cellType %>% filter(condition == "scrambled") %>% pull(auprc)

wilcox.test(auprc_classifier, auprc_scrambledClassifier, paired = TRUE) #V = 5050, p-value < 2.2e-16 #20240517 MEM
t.test(auprc_classifier, auprc_scrambledClassifier, paired = TRUE) #t = 48.79, df = 99, p-value < 2.2e-16, #20240517 MEM

# AUROC
auroc_classifier <- scrambledClassifier_cellType %>% filter(condition == "classifier") %>% pull(auroc)
auroc_scrambledClassifier <- scrambledClassifier_cellType %>% filter(condition == "scrambled") %>% pull(auroc)

wilcox.test(auroc_classifier, auroc_scrambledClassifier, paired = TRUE) #V = 5050, p-value < 2.2e-16 #20240517 MEM
t.test(auroc_classifier, auroc_scrambledClassifier, paired = TRUE) #t = 76.525, df = 99, p-value < 2.2e-16, #20240517 MEM




# SCRAMBLED TEST SAMPLE

dataset = "ClonMapper"
sample1_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample1Scrambled_sample2Classifier.csv"))
sample2_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample2Scrambled_sample1Classifier.csv"))
clonMapper_scrambTest = bind_rows(sample1_scrambTest, sample2_scrambTest)

dataset = "TREX"
sample1_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample1Scrambled_sample2Classifier.csv"))
sample2_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample2Scrambled_sample1Classifier.csv"))
trex_scrambTest = bind_rows(sample1_scrambTest, sample2_scrambTest)

dataset = "SPLINTR"
sample1_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample1Scrambled_sample2Classifier.csv"))
sample2_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample2Scrambled_sample1Classifier.csv"))
sample3_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample3Scrambled_sample4Classifier.csv"))
sample4_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample4Scrambled_sample3Classifier.csv"))
splintr_scrambTest = bind_rows(sample1_scrambTest, sample2_scrambTest, sample3_scrambTest, sample4_scrambTest)

dataset = "FateMap" #Goyal et al. 1
sample1_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample1Scrambled_sample2Classifier.csv"))
sample2_scrambTest = read.csv(paste0(dataset_dir, dataset, "/sample2Scrambled_sample1Classifier.csv"))
fateMap_scrambTest = bind_rows(sample1_scrambTest, sample2_scrambTest)

scrambled_cellType = bind_rows(clonMapper_scrambTest, trex_scrambTest, splintr_scrambTest, fateMap_scrambTest) %>% mutate(scrambled = TRUE)

allCellTypeClassifiers_forScrambled = cellTypeSpecificClassifiers %>% filter(classifier == TRUE) %>% mutate(scrambled = FALSE)
scrambled_cellType = scrambled_cellType %>% bind_rows(allCellTypeClassifiers_forScrambled) %>% select(-X, -best_params, -row_number, -sample, -dbl_exp, -dbl_act, -technology_dataset, -technology, -classifier)

#write.csv(scrambled_cellType, file = paste0(dataset_dir, "allCellTypeClassifiersScrambledResults.csv")) #20240517 MEM

# AUPRC
plot = ggplot(scrambled_cellType, aes(x = scrambled, y = auprc, color = scrambled)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUPRC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auprc_scrambled.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auprc_scrambled.png'), width = 2, height = 4) #20240517 MEM

# AUROC
plot = ggplot(scrambled_cellType, aes(x = scrambled, y = auroc, color = scrambled)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUROC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auroc_scrambled.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allCellTypeClassifiers_auroc_scrambled.png'), width = 2, height = 4) #20240517 MEM


### STATISTICS

# AUPRC
auprc_notScrambled <- scrambled_cellType %>% filter(scrambled == FALSE) %>% pull(auprc)
auprc_scrambled <- scrambled_cellType %>% filter(scrambled == TRUE) %>% pull(auprc)

wilcox.test(auprc_notScrambled, auprc_scrambled, paired = TRUE) #V = 5050, p-value < 2.2e-16 #20240517 MEM
t.test(auprc_notScrambled, auprc_scrambled, paired = TRUE) #t = 32.46, df = 99, p-value < 2.2e-16, #20240517 MEM

# AUROC
auroc_notScrambled <- scrambled_cellType %>% filter(scrambled == FALSE) %>% pull(auroc)
auroc_scrambled <- scrambled_cellType %>% filter(scrambled == TRUE) %>% pull(auroc)

wilcox.test(auroc_notScrambled, auroc_scrambled, paired = TRUE) #V = 5050, p-value < 2.2e-16 #20240517 MEM
t.test(auroc_notScrambled, auroc_scrambled, paired = TRUE) #t = 69.76, df = 99, p-value < 2.2e-16, #20240517 MEM










### for mouse and human integrated data ###

# SCRAMBLED CLASSIFIER 

A_human_scrambClass = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/A_human.csv")
B_human_scrambClass = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/B_human.csv")
A_mouse_scrambClass = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/A_mouse.csv")
B_mouse_scrambClass = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/B_mouse.csv")

scrambledClassifier_HM = bind_rows(A_human_scrambClass, B_human_scrambClass, A_mouse_scrambClass, B_mouse_scrambClass) %>% mutate(condition = str_replace(condition, ".*neg_control_scrambled.*", "scrambled"),
                                                                                                                                  condition = if_else(condition == as.character(dataset), "classifier", condition))
scrambledClassifier_HM = scrambledClassifier_HM  %>% filter(condition %in% c("scrambled", "classifier"))

#write.csv(scrambledClassifier_HM, file = paste0(dataset_dir, "allHumanMouseClassifiersScrambledClassifierResults.csv")) #20240517 MEM

# AUPRC
plot = ggplot(scrambledClassifier_HM, aes(x = condition, y = auprc, color = condition)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUPRC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auprc_scrambledClassifier.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auprc_scrambledClassifier.png'), width = 2, height = 4) #20240517 MEM

# AUROC
plot = ggplot(scrambledClassifier_HM, aes(x = condition, y = auroc, color = condition)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUROC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auroc_scrambledClassifier.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auroc_scrambledClassifier.png'), width = 2, height = 4) #20240517 MEM

### STATISTICS

# AUPRC
auprc_classifier <- scrambledClassifier_HM %>% filter(condition == "classifier") %>% pull(auprc)
auprc_scrambledClassifier <- scrambledClassifier_HM %>% filter(condition == "scrambled") %>% pull(auprc)

wilcox.test(auprc_classifier, auprc_scrambledClassifier, paired = TRUE) #V = 820, p-value = 1.819e-12 #20240517 MEM
t.test(auprc_classifier, auprc_scrambledClassifier, paired = TRUE) #t = 36.31, df = 39, p-value < 2.2e-16, #20240517 MEM

# AUROC
auroc_classifier <- scrambledClassifier_HM %>% filter(condition == "classifier") %>% pull(auroc)
auroc_scrambledClassifier <- scrambledClassifier_HM %>% filter(condition == "scrambled") %>% pull(auroc)

wilcox.test(auroc_classifier, auroc_scrambledClassifier, paired = TRUE) #V = 820, p-value = 1.819e-12 #20240517 MEM
t.test(auroc_classifier, auroc_scrambledClassifier, paired = TRUE) #t = 71.727, df = 39, p-value < 2.2e-16, #20240517 MEM



# SCRAMBLED TEST SAMPLE

A_human_FM2 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/FateMap2_B_humanScrambled_A_humanClassifier.csv")
A_human_CM2 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/ClonMapper2_B_humanScrambled_A_humanClassifier.csv")
B_human_FM1 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/FateMap1_A_humanScrambled_B_humanClassifier.csv")
B_human_CM1 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/ClonMapper1_A_humanScrambled_B_humanClassifier.csv")

scrambled_human = bind_rows(A_human_FM2, A_human_CM2, B_human_FM1, B_human_CM1) %>% mutate(technology = "human_all")

A_mouse_SP2 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/SPLINTR2_B_mouseScrambled_A_mouseClassifier.csv")
A_mouse_SP4 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/SPLINTR4_B_mouseScrambled_A_mouseClassifier.csv")
A_mouse_TX2 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/TREX2_B_mouseScrambled_A_mouseClassifier.csv")
B_mouse_SP1 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/SPLINTR1_A_mouseScrambled_B_mouseClassifier.csv")
B_mouse_SP3 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/SPLINTR3_A_mouseScrambled_B_mouseClassifier.csv")
B_mouse_TX1 = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/classifier/crossExperiment/results/datasetSummaries/all/TREX1_A_mouseScrambled_B_mouseClassifier.csv")

scrambled_mouse = bind_rows(A_mouse_SP2, A_mouse_SP4, A_mouse_TX2, B_mouse_SP1, B_mouse_SP3, B_mouse_TX1) %>% mutate(technology = "mouse_all")

scrambled_HM = bind_rows(scrambled_human, scrambled_mouse) %>% mutate(scrambled = TRUE)
allHumanMouseClassifiers_forScrambled = allHumanMouseClassifiers %>% filter(classifier == TRUE) %>% mutate(scrambled = FALSE)
scrambled_HM = scrambled_HM %>% bind_rows(allHumanMouseClassifiers_forScrambled) %>% select(-best_params, -classifier, -sample, -dbl_exp, -dbl_act, -technology_dataset)

#write.csv(scrambled_HM, file = paste0(dataset_dir, "allHumanMouseClassifiersScrambledResults.csv")) #20240516 MEM

# AUPRC
plot = ggplot(scrambled_HM, aes(x = scrambled, y = auprc, color = scrambled)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUPRC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auprc_scrambled.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auprc_scrambled.png'), width = 2, height = 4) #20240517 MEM

# AUROC
plot = ggplot(scrambled_HM, aes(x = scrambled, y = auroc, color = scrambled)) +
  #geom_pointrange(data = allClassifiers_summary, aes(x = technology, y =avgAuprc, ymin = avgAuprc-sdAuprc, ymax = avgAuprc+sdAuprc, fill = "condition"), position = position_dodge(0.6), color = "black", size = 1, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5, shape = 16) +  
  #scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0,1.01)) +
  labs( y = "AUROC", x = "") +
  theme_classic() +
  theme(legend.position = "none")
plot
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auroc_scrambled.svg'), width = 2, height = 4) #20240517 MEM
#ggsave(plot, file = paste0(plot_directory, 'allHumanMouseClassifiers_auroc_scrambled.png'), width = 2, height = 4) #20240517 MEM


### STATISTICS

# AUPRC
auprc_notScrambled <- scrambled_HM %>% filter(scrambled == FALSE) %>% pull(auprc)
auprc_scrambled <- scrambled_HM %>% filter(scrambled == TRUE) %>% pull(auprc)

wilcox.test(auprc_notScrambled, auprc_scrambled, paired = TRUE) #V = 5050, p-value < 2.2e-16 #20240517 MEM
t.test(auprc_notScrambled, auprc_scrambled, paired = TRUE) #t = 37.373, df = 99, p-value < 2.2e-16, #20240517 MEM

# AUROC
auroc_notScrambled <- scrambled_HM %>% filter(scrambled == FALSE) %>% pull(auroc)
auroc_scrambled <- scrambled_HM %>% filter(scrambled == TRUE) %>% pull(auroc)

wilcox.test(auroc_notScrambled, auroc_scrambled, paired = TRUE) #V = 5050, p-value < 2.2e-16 #20240517 MEM
t.test(auroc_notScrambled, auroc_scrambled, paired = TRUE) #t = 89.652, df = 99, p-value < 2.2e-16, #20240517 MEM








