### plotting differential expression results together for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20240206
### Last edited by Madeline E Melzer on 20240215

library(tidyverse)
library(ggplot2)
library(dplyr)
library(glue)
library(readr)
library(PRROC)
library(conflicted)
conflict_prefer("filter", "dplyr")

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


plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/differentialExpression/"
deData = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/DE/plotData/de.csv") #this is from Charles

################################
#                              #
#           Wilcox             #
#                              #
################################


wilcox <- deData %>% filter(method == "wilcox") %>% mutate()
wilcox$dataset = as_factor(wilcox$dataset)

wilcox = wilcox %>% dplyr::filter(dataset %in% barcoded)
wilcox$dataset <- factor(wilcox$dataset, levels = barcoded)
wilcox <- wilcox %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

plot = ggplot(wilcox, aes(x = dataset, y = precision)) +
  geom_line(data = wilcox, aes(x = dataset, y = precision, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'DE_Wilcox_precision.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'DE_Wilcox_precision.png'), width = 4, height = 3)

plot = ggplot(wilcox, aes(x = dataset, y = recall)) +
  geom_line(data = wilcox, aes(x = dataset, y = recall, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'DE_Wilcox_recall.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'DE_Wilcox_recall.png'), width = 4, height = 3)

plot = ggplot(wilcox, aes(x = dataset, y = TNR)) +
  geom_line(data = wilcox, aes(x = dataset, y = TNR, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'DE_Wilcox_TNR.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'DE_Wilcox_TNR.png'), width = 4, height = 3)




### wilcox statistics 20240321

wilcox_precision = wilcox %>%
  pivot_wider(names_from = dbl_pct, values_from = precision, id_cols = c(method, dataset))
 
wilcox.test(wilcox_precision$"10%", wilcox_precision$"20%", paired = TRUE) # V = 120, p-value = 6.104e-05
wilcox.test(wilcox_precision$"10%", wilcox_precision$"40%", paired = TRUE) #V = 120, p-value = 6.104e-05

wilcox_recall = wilcox %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(method, dataset))

wilcox.test(wilcox_recall$"10%", wilcox_recall$"20%", paired = TRUE) #V = 105, p-value = 0.001097
wilcox.test(wilcox_recall$"10%", wilcox_recall$"40%", paired = TRUE, exact = TRUE) #V = 120, p-value = 6.104e-05

wilcox_TNR = wilcox %>%
  pivot_wider(names_from = dbl_pct, values_from = TNR, id_cols = c(method, dataset))

wilcox.test(wilcox_TNR$"10%", wilcox_TNR$"20%", paired = TRUE, exact = TRUE) #V = 120, p-value = 6.104e-05
wilcox.test(wilcox_TNR$"10%", wilcox_TNR$"40%", paired = TRUE, exact = TRUE) #V = 120, p-value = 6.104e-05


# t test

wilcox_precision = wilcox %>%
  pivot_wider(names_from = dbl_pct, values_from = precision, id_cols = c(method, dataset))

t.test(wilcox_precision$"10%", wilcox_precision$"20%", paired = TRUE) # t = 7.1597, df = 14, p-value = 4.861e-06
t.test(wilcox_precision$"20%", wilcox_precision$"40%", paired = TRUE) # t = 11.897, df = 14, p-value = 1.044e-08

wilcox_recall = wilcox %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(method, dataset))

t.test(wilcox_recall$"10%", wilcox_recall$"20%", paired = TRUE) # t = 4.6624, df = 14, p-value = 0.0003663
t.test(wilcox_recall$"20%", wilcox_recall$"40%", paired = TRUE) # t = 6.5495, df = 14, p-value = 1.292e-05

wilcox_TNR = wilcox %>%
  pivot_wider(names_from = dbl_pct, values_from = TNR, id_cols = c(method, dataset))

t.test(wilcox_TNR$"10%", wilcox_TNR$"20%", paired = TRUE) # t = 5.5675, df = 14, p-value = 6.935e-05
t.test(wilcox_TNR$"20%", wilcox_TNR$"40%", paired = TRUE) # t = 5.015, df = 14, p-value = 0.0001892







################################
#                              #
#             MAST             #
#                              #
################################

mast <- deData %>% filter(method == "MAST")
mast$dataset = as_factor(mast$dataset)

mast = mast %>% dplyr::filter(dataset %in% barcoded)
mast$dataset <- factor(mast$dataset, levels = barcoded)
mast <- mast %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

plot = ggplot(mast, aes(x = dataset, y = precision)) +
  geom_line(data = mast, aes(x = dataset, y = precision, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'DE_MAST_precision.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'DE_MAST_precision.png'), width = 4, height = 3)

plot = ggplot(mast, aes(x = dataset, y = recall)) +
  geom_line(data = mast, aes(x = dataset, y = recall, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'DE_MAST_recall.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'DE_MAST_recall.png'), width = 4, height = 3)

plot = ggplot(mast, aes(x = dataset, y = TNR)) +
  geom_line(data = mast, aes(x = dataset, y = TNR, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'DE_MAST_TNR.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'DE_MAST_TNR.png'), width = 4, height = 3)





### MAST statistics 20240321

mast_precision = mast %>%
  pivot_wider(names_from = dbl_pct, values_from = precision, id_cols = c(method, dataset))

wilcox.test(mast_precision$"10%", mast_precision$"20%", paired = TRUE, exact = TRUE) # V = 120, p-value = 6.104e-05
wilcox.test(mast_precision$"10%", mast_precision$"40%", paired = TRUE) #V = 120, p-value = 6.104e-05


mast_recall = mast %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(method, dataset))

wilcox.test(mast_recall$"10%", mast_recall$"20%", paired = TRUE) #V = 120, p-value = 6.104e-05
wilcox.test(mast_recall$"10%", mast_recall$"40%", paired = TRUE) #V = 120, p-value = 6.104e-05

mast_TNR = mast %>%
  pivot_wider(names_from = dbl_pct, values_from = TNR, id_cols = c(method, dataset))

wilcox.test(mast_TNR$"10%", mast_TNR$"20%", paired = TRUE) #V = 120, p-value = 6.104e-05
wilcox.test(mast_TNR$"10%", mast_TNR$"40%", paired = TRUE) #V = 120, p-value = 6.104e-05



# t test, 20240321

t.test(mast_precision$"10%", mast_precision$"20%", paired = TRUE) # t = 5.7094, df = 14, p-value = 5.396e-05
t.test(mast_precision$"20%", mast_precision$"40%", paired = TRUE) #t = 10.619, df = 14, p-value = 4.419e-08

mast_recall = mast %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(method, dataset))

t.test(mast_recall$"10%", mast_recall$"20%", paired = TRUE) #t = 5.2538, df = 14, p-value = 0.000122
t.test(mast_recall$"20%", mast_recall$"40%", paired = TRUE) #t = 7.2112, df = 14, p-value = 4.486e-06

mast_TNR = mast %>%
  pivot_wider(names_from = dbl_pct, values_from = TNR, id_cols = c(method, dataset))

t.test(mast_TNR$"10%", mast_TNR$"20%", paired = TRUE) #t = 5.1003, df = 14, p-value = 0.0001616
t.test(mast_TNR$"20%", mast_TNR$"40%", paired = TRUE) #t = 5.1429, df = 14, p-value = 0.0001494













#all of this is useless !!!! I did it before I saw that Charles has alredy done the pr and roc calculations
subfolders <- list.dirs(dataDirectory, recursive = FALSE)

deData <- map_dfr(subfolders, ~{
  detection_files <- list.files(.x, pattern = "*.csv", full.names = TRUE)
  folder_name = basename(.x)
  # Iterate over each file and read + process the data
  map_dfr(detection_files, ~{
    file_name <- basename(.x)
    parts <- str_split(file_name, "_|\\.csv", simplify = TRUE)
    doubletRate_raw = parts[1]
    
    read_csv(.x) %>%
      mutate(dataset = folder_name,  # Use dirname() to get the folder name
             clusterA = str_extract(file_name, "(?<=_)[0-9]+(?=_vs\\.)"),  # Number after the first underscore and before '_vs.'
             clusterB = str_extract(file_name, "(?<=_vs\\._)[0-9]+"),  # Number after '_vs._'
             method = str_extract(file_name, "wilcox|MAST"),
             doubletRate = case_when(
               str_detect(doubletRate_raw, "control") ~ "0",
               str_detect(doubletRate_raw, "twenty") ~ "20",
               str_detect(doubletRate_raw, "forty") ~ "40",
               str_detect(doubletRate_raw, "ten") ~ "10",
               TRUE ~ doubletRate_raw # Handles any other unexpected values
             )
      )
  })
})

write_csv(deData, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/DE/differentialExpression.csv")

de_rocpr <- tibble(dataset = character(), method = character(), doubletRate = character(), roc = numeric(), pr = numeric())

# Iterate over combinations of dataset, method, and doubletRate
for (dataset in unique(deData$dataset)) {
  print(dataset)
  for (method in unique(deData$method)) {
    print(method)
    control_scores <- deData %>%
      filter(dataset == dataset, method == method, doubletRate == "0") %>%
      pull(avg_log2FC)
    
    for (rate in unique(deData$doubletRate)) {
      print(rate)
      if (rate != "0") {
        exp_scores <- deData %>%
          filter(dataset == dataset, method == method, doubletRate == rate) %>%
          pull(avg_log2FC)
        
        print(length(control_scores))
        print(length(exp_scores))
        
        pr <- pr.curve(scores.class0 = control_scores, scores.class1 = exp_scores, curve = TRUE)
        cur_pr <- pr$auc.integral
        roc <- roc.curve(scores.class0 = control_scores, scores.class1 = exp_scores, curve = TRUE)
        cur_roc <- roc$auc
        
        # Append results to the dataframe
        de_rocpr <- de_rocpr %>%
          add_row(dataset = dataset, method = method, doubletRate = rate, 
                  roc = cur_roc, pr = cur_pr)
      }
    }
  }
}
write_csv(de_rocpr, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/DE/differentialExpression_rocpr.csv")
colnames(de_rocpr)[colnames(de_rocpr) == "roc"] <- "auroc"
colnames(de_rocpr)[colnames(de_rocpr) == "pr"] <- "auprc"
