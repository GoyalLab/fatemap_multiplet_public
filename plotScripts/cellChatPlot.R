### plotting cell-cell communication- Cell Chat results together for Zhang Melzer et al 2024
### Created by Madeline E Melzer on 20240215
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


plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/plots/cellChat/"
ccData = read_csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/cellChat/cc.csv") #this is from Charles

################################
#                              #
#          Cell Chat           #
#                              #
################################

ccData = ccData %>% dplyr::filter(dataset %in% barcoded)
ccData$dataset <- factor(ccData$dataset, levels = barcoded)
ccData <- ccData %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

write.csv(ccData, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/cellChat/cc_formatted.csv")


plot = ggplot(ccData, aes(x = dataset, y = precision)) +
  geom_line(data = ccData, aes(x = dataset, y = precision, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'cellChat_precision.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'cellChat_precision.png'), width = 4, height = 3)

plot = ggplot(ccData, aes(x = dataset, y = recall)) +
  geom_line(data = ccData, aes(x = dataset, y = recall, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot

ggsave(plot, file = paste0(plotDirectory, 'cellChat_recall.svg'), width = 4, height = 3)
ggsave(plot, file = paste0(plotDirectory, 'cellChat_recall.png'), width = 4, height = 3)




### cellChat statistics 20240321

ccData_precision = ccData %>%
  pivot_wider(names_from = dbl_pct, values_from = precision, id_cols = c(dataset))

wilcox.test(ccData_precision$"10%", ccData_precision$"20%", paired = TRUE) # V = 120, p-value = 6.104e-05
wilcox.test(ccData_precision$"10%", ccData_precision$"40%", paired = TRUE) #V = 120, p-value = 6.104e-05

ccData_recall = ccData %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(dataset))

wilcox.test(ccData_recall$"10%", ccData_recall$"20%", paired = TRUE) #V = 99, p-value = 0.00388
wilcox.test(ccData_recall$"10%", ccData_recall$"40%", paired = TRUE) #V = 120, p-value = 6.104e-05


# t test

ccData_precision = ccData %>%
  pivot_wider(names_from = dbl_pct, values_from = precision, id_cols = c(dataset))

t.test(ccData_precision$"10%", ccData_precision$"20%", paired = TRUE) # t = 8.6363, df = 14, p-value = 5.557e-07
t.test(ccData_precision$"20%", ccData_precision$"40%", paired = TRUE) # t = 8.2019, df = 14, p-value = 1.023e-06

ccData_recall = ccData %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(dataset))

t.test(ccData_recall$"10%", ccData_recall$"20%", paired = TRUE) # t = 3.9325, df = 14, p-value = 0.001503
t.test(ccData_recall$"20%", ccData_recall$"40%", paired = TRUE) # t = 4.3833, df = 14, p-value = 0.0006246


################################
#                              #
#          CellPhoneDB        #
#          #20240516                    #
################################

cpdbData = read.csv("/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/cellChat/cpdb.csv")


cpdbData = cpdbData %>% dplyr::filter(dataset %in% barcoded)
cpdbData$dataset <- factor(cpdbData$dataset, levels = barcoded)
cpdbData <- cpdbData %>%
  mutate(dataset = recode(dataset, !!!barcodedRename))

write.csv(cpdbData, "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/charles/cellChat/cpdb_formatted.csv")

plot = ggplot(cpdbData, aes(x = dataset, y = precision)) +
  geom_line(data = cpdbData, aes(x = dataset, y = precision, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plotDirectory, 'cellPhoneDB_precision.svg'), width = 4, height = 3) #20240516 MEM
#ggsave(plot, file = paste0(plotDirectory, 'cellPhoneDB_precision.png'), width = 4, height = 3) #20240516 MEM

plot = ggplot(cpdbData, aes(x = dataset, y = recall)) +
  geom_line(data = cpdbData, aes(x = dataset, y = recall, group = dbl_pct, color = as.factor(dbl_pct)), size = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1.01)) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot
#ggsave(plot, file = paste0(plotDirectory, 'cellPhoneDB_recall.svg'), width = 4, height = 3) #20240516 MEM
#ggsave(plot, file = paste0(plotDirectory, 'cellPhoneDB_recall.png'), width = 4, height = 3) #20240516 MEM




### cellPhoneDB statistics 20240516 MEM

cpdbData_precision = cpdbData %>%
  pivot_wider(names_from = dbl_pct, values_from = precision, id_cols = c(dataset))

wilcox.test(cpdbData_precision$"10%", cpdbData_precision$"20%", paired = TRUE) # V = 66, p-value = 0.0009766
wilcox.test(cpdbData_precision$"10%", cpdbData_precision$"40%", paired = TRUE) #V = 66, p-value = 0.0009766

cpdbData_recall = cpdbData %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(dataset))

wilcox.test(cpdbData_recall$"10%", cpdbData_recall$"20%", paired = TRUE) #V = 25, p-value = 0.8385
wilcox.test(cpdbData_recall$"10%", cpdbData_recall$"40%", paired = TRUE) #V = 29, p-value = 0.9188


# t test

cpdbData_precision = cpdbData %>%
  pivot_wider(names_from = dbl_pct, values_from = precision, id_cols = c(dataset))

t.test(cpdbData_precision$"10%", cpdbData_precision$"20%", paired = TRUE) # t = 6.8282, df = 10, p-value = 4.579e-05
t.test(cpdbData_precision$"20%", cpdbData_precision$"40%", paired = TRUE) # t = 14.542, df = 10, p-value = 4.707e-08

cpdbData_recall = cpdbData %>%
  pivot_wider(names_from = dbl_pct, values_from = recall, id_cols = c(dataset))

t.test(cpdbData_recall$"10%", cpdbData_recall$"20%", paired = TRUE) # t = -0.26601, df = 10, p-value = 0.7956
t.test(cpdbData_recall$"20%", cpdbData_recall$"40%", paired = TRUE) # t = 0.65048, df = 10, p-value = 0.53




