library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library("ggsci")
library(forcats)
library(stringr)

dataDirectory <- "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plotData/heterogenityAnalysis/"

plotDirectory = "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plotData/heterogenityAnalysis/plots/"

HetAll <- tibble(read.table(file = paste0(dataDirectory,"all_average.csv"), header = TRUE, sep = ",", stringsAsFactors=F)) %>% select (-diff)

HetAllLong = HetAll %>% gather("condition", "AUPRC", -dataset_id, -method) %>% mutate(pointClassifier = paste0(method, dataset_id))

HetAllLongSummary = HetAllLong %>% group_by(condition, method) %>% summarise(meanAUPRC = mean(AUPRC)) 

plot1 = ggplot(data=HetAllLong, aes(x=condition, y=AUPRC)) + 
  geom_point(aes(color = method), alpha=0.2) + geom_line(aes(group = pointClassifier), alpha=0.2) +
  geom_point(data = HetAllLongSummary, aes(x=condition, y=meanAUPRC, color = method), size = 3) +
  geom_line(data = HetAllLongSummary, aes(x=condition, y=meanAUPRC, group = method)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~method)

HetAllDoubletCell = HetAll %>% filter(method == "doublet_cell")  
HetAllDoubletFinder = HetAll %>% filter(method == "doublet_finder")  

wilcox.test(HetAll$adj_test, HetAll$all_cluster, paired = TRUE)

wilcox.test(HetAllDoubletCell$adj_test, HetAllDoubletCell$all_cluster, paired = TRUE)

wilcox.test(HetAllDoubletFinder$adj_test, HetAllDoubletFinder$all_cluster, paired = TRUE)

#######All datasets Together#######

HetNonAverage <- tibble(read.table(file = paste0(dataDirectory,"all.tsv"), header = TRUE, stringsAsFactors=F)) %>% filter(type != "three_cluster")

HetNonAverage = HetNonAverage %>% group_by(method, type, dataset_id) %>% mutate(replicateID = c("1", "2", "3")) %>% mutate(pointClassifier = paste0(method, dataset_id, replicateID))

HetNonAverageSummary = HetNonAverage %>% group_by(type, method) %>% summarise(meanPR = mean(pr),
                                                                                meanROC = mean(roc)) 

plot1 = ggplot(data=HetNonAverage, aes(x=type, y=pr)) + 
  geom_point(aes(color = method), alpha=0.2) + geom_line(aes(group = pointClassifier), alpha=0.2) +
  geom_point(data = HetNonAverageSummary, aes(x=type, y=meanPR, color = method), size = 3) +
  geom_line(data = HetNonAverageSummary, aes(x=type, y=meanPR, group = method)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~method)

HetNonAverageSpread = HetNonAverage %>% ungroup() %>% select(pr,pointClassifier,type, method) 
HetNonAverageSpreadDoublet_cell = HetNonAverageSpread %>% spread(type,pr) %>% filter(method =="doublet_cell")
HetNonAverageSpreadDoublet_finder = HetNonAverageSpread %>% spread(type,pr) %>% filter(method =="doublet_finder")


wilcox.test(HetNonAverageSpreadDoublet_cell$adj_test, HetNonAverageSpreadDoublet_cell$all_cluster, paired = TRUE)
wilcox.test(HetNonAverageSpreadDoublet_finder$adj_test, HetNonAverageSpreadDoublet_finder$all_cluster, paired = TRUE)

HetNonAverageSpread = HetNonAverage %>% ungroup() %>% select(roc,pointClassifier,type, method) 
HetNonAverageSpreadDoublet_cell = HetNonAverageSpread %>% spread(type,roc) %>% filter(method =="doublet_cell")
HetNonAverageSpreadDoublet_finder = HetNonAverageSpread %>% spread(type,roc) %>% filter(method =="doublet_finder")

plot1 = ggplot(data=HetNonAverage, aes(x=type, y=roc)) + 
  geom_point(aes(color = method), alpha=0.2) + geom_line(aes(group = pointClassifier), alpha=0.2) +
  geom_point(data = HetNonAverageSummary, aes(x=type, y=meanROC, color = method), size = 3) +
  geom_line(data = HetNonAverageSummary, aes(x=type, y=meanROC, group = method)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~method)

######################################################################################################
################FinalAnalysis for the paper: heterogenity adjacent vs non ############################
######################################################################################################

dataDirectory <- "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plotData/heterogenityAnalysis/"
plotDirectory = "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plots/heterogenityAnalysis/"

HetNonAverage <- tibble(read.table(file = paste0(dataDirectory,"all.tsv"), header = TRUE, stringsAsFactors=F)) %>% filter(type != "three_cluster")

HetNonAverage = HetNonAverage %>% group_by(method, type, dataset_id) %>% summarise(pr = mean(pr),roc = mean(roc))
HetNonAverage = HetNonAverage %>% mutate(pointClassifier = paste0(method, dataset_id)) 

HetNonAverageSummary = HetNonAverage %>% group_by(type, method) %>% summarise(meanPR = mean(pr),
                                                                              meanROC = mean(roc)) 


plotPR = ggplot(data=HetNonAverage, aes(x=type, y=pr)) + 
  geom_point(aes(color = method), alpha=0.2) + geom_line(aes(group = pointClassifier), alpha=0.2) +
  geom_point(data = HetNonAverageSummary, aes(x=type, y=meanPR, color = method), size = 3) +
  geom_line(data = HetNonAverageSummary, aes(x=type, y=meanPR, group = method)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~method)

ggsave(plotPR, file = paste0(plotDirectory, 'individualHeterogenity_PR.svg'), width = 5, height = 4)


plotROC = ggplot(data=HetNonAverage, aes(x=type, y=roc)) + 
  geom_point(aes(color = method), alpha=0.2) + geom_line(aes(group = pointClassifier), alpha=0.2) +
  geom_point(data = HetNonAverageSummary, aes(x=type, y=meanROC, color = method), size = 3) +
  geom_line(data = HetNonAverageSummary, aes(x=type, y=meanROC, group = method)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~method)

ggsave(plotROC, file = paste0(plotDirectory, 'individualHeterogenity_ROC.svg'), width = 5, height = 4)


HetNonAverageSpread = HetNonAverage %>% ungroup() %>% select(pr,pointClassifier,type, method) 
HetNonAverageSpreadDoublet_cell = HetNonAverageSpread %>% spread(type,pr) %>% filter(method =="doublet_cell")
HetNonAverageSpreadDoublet_finder = HetNonAverageSpread %>% spread(type,pr) %>% filter(method =="doublet_finder")


wilcox.test(HetNonAverageSpreadDoublet_cell$adj_test, HetNonAverageSpreadDoublet_cell$all_cluster, paired = TRUE)
wilcox.test(HetNonAverageSpreadDoublet_finder$adj_test, HetNonAverageSpreadDoublet_finder$all_cluster, paired = TRUE)



HetNonAverageSpread = HetNonAverage %>% ungroup() %>% select(roc,pointClassifier,type, method) 
HetNonAverageSpreadDoublet_cell = HetNonAverageSpread %>% spread(type,roc) %>% filter(method =="doublet_cell")
HetNonAverageSpreadDoublet_finder = HetNonAverageSpread %>% spread(type,roc) %>% filter(method =="doublet_finder")


wilcox.test(HetNonAverageSpreadDoublet_cell$adj_test, HetNonAverageSpreadDoublet_cell$all_cluster, paired = TRUE)
wilcox.test(HetNonAverageSpreadDoublet_finder$adj_test, HetNonAverageSpreadDoublet_finder$all_cluster, paired = TRUE)

################################################################################################
####################################ranking analysis cxds######################################
################################################################################################

dataDirectory <- "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plotData/similarityScore/"

plotDirectory = "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plots/similarityScore/"

data = tibble(read.table(file = paste0(dataDirectory,"allPoints_similarityScores.csv"), header = TRUE, sep = ",", stringsAsFactors=F))
data1 = tibble(read.table(file = paste0(dataDirectory,"allPoints_similarityScores.csv"), header = TRUE, sep = ",", stringsAsFactors=F))

data <- data %>% 
  mutate(dataType = ifelse(methodA == "cxds" | methodB == "cxds", "cxds", "other"))
  
scale_rank <- function(df) {
  df$rank <- rank(df$value, ties.method= "first")  # Rank by area (descending)
  max_rank <- max(df$rank)
  df$scaled_rank <- (df$rank - 1) / (max_rank - 1) * 100
  return(df)
}

dataRank <- scale_rank(data)

nrow(dataRank %>% filter(dataType == "cxds"))
nrow(dataRank %>% filter(dataType == "other"))

nrow(dataRank %>% filter(rank<201) %>% filter(dataType == "cxds"))
nrow(dataRank %>% filter(rank<201) %>% filter(dataType == "other"))

plot = ggplot(dataRank, aes(x = rank, y = value, fill = as.factor(dataType))) +
  geom_col() +
  geom_rug(sides="b", mapping = aes(color = as.factor(dataType))) +
  scale_color_manual(values = c("hotpink3", "lightgray")) +
  scale_fill_manual(values = c("hotpink3", "lightgray")) +
  theme_classic() +
  theme(legend.position = "none")
ggsave(plot, file = paste0(plotDirectory, 'cxds_Ranking_v1.svg'), width = 6, height = 4)


################################################################################################
####################################TNR for Amulet ######################################
################################################################################################

dataDirectory <- "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plotData/AmuletAnalysis/"
plotDirectory = "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plots/AmuletAnalysis/"

dataTNR = tibble(read.table(file = paste0(dataDirectory,"combined_umi1_res.csv"), header = TRUE, sep = ",", stringsAsFactors=F))

TNRAverage = dataTNR %>% summarise(meanTNR = mean(TNR),
                                   sdTNR = sd(TNR))


plot = ggplot() +
  geom_jitter(data = dataTNR, aes(x="age", y=TNR),size=2, shape = 16, width = 0.2, height = 0.01) +
  geom_pointrange(data = TNRAverage, aes(x = "age", y =meanTNR, ymin = meanTNR-sdTNR, ymax = meanTNR+sdTNR)) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  ylim(0,1)

ggsave(plot, file = paste0(plotDirectory, 'TNR_Amulet.svg'), width = 2, height = 5)

################################################################################################
####################################average vs summed ##########################################
################################################################################################

dataDirectory <- "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plotData/averagedAndSummedDoublets/"
plotDirectory = "/Volumes/GoogleDrive/My Drive/Northwestern/papers/ZhangMelzerEtAl/Revisions/plots/averagedAndSummedDoublets/"

dataSumAverage = tibble(read.table(file = paste0(dataDirectory,"averagedAndSummedDoublets.csv"), header = TRUE, sep = ",", stringsAsFactors=F))

dataSumAverage <- dataSunAverage %>% mutate(pointClassifier = paste0(dataset, condition, sample))

dataSunAverage.summary <- dataSunAverage %>% group_by(condition, type) %>% summarise(auprc_mean = mean(auprc),
                                                                                     auroc_mean = mean(auroc),
                                                                                     TNRmean = mean(TNR))

dataSunAverage.summary.summary <- dataSunAverage.summary %>% group_by(type) %>% summarise(auprc_mean = mean(auprc_mean),
                                                                                  auroc_mean = mean(auroc_mean),
                                                                                  TNRmean = mean(TNRmean)) 

plotPR = ggplot(data=dataSumAverage, aes(x=type, y=auprc)) + 
  geom_line(aes(color = condition, group = pointClassifier), alpha=0.1) +
  geom_point(data = dataSunAverage.summary, aes(x=type, y=auprc_mean), size = 3) +
  geom_line(data = dataSunAverage.summary, aes(x=type, y=auprc_mean, group = condition)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~condition, ncol = 4)

ggsave(plotPR, file = paste0(plotDirectory, 'AverageSum_PR.svg'), width = 6, height = 5)

plotROC = ggplot(data=dataSumAverage, aes(x=type, y=auroc)) + 
  geom_line(aes(color = condition, group = pointClassifier), alpha=0.1) +
  geom_point(data = dataSunAverage.summary, aes(x=type, y=auroc_mean), size = 3) +
  geom_line(data = dataSunAverage.summary, aes(x=type, y=auroc_mean, group = condition)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~condition, ncol = 4)

ggsave(plotROC, file = paste0(plotDirectory, 'AverageSum_ROC.svg'), width = 6, height = 5)

plotTNR = ggplot(data=dataSumAverage, aes(x=type, y=TNR)) + 
  geom_line(aes(color = condition, group = pointClassifier), alpha=0.1) +
  geom_point(data = dataSunAverage.summary, aes(x=type, y=TNRmean), size = 3) +
  geom_line(data = dataSunAverage.summary, aes(x=type, y=TNRmean, group = condition)) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  facet_wrap(~condition, ncol = 4)

ggsave(plotTNR, file = paste0(plotDirectory, 'AverageSum_TNR.svg'), width = 6, height = 5)

################################################################################################
#################################### Stat testing for average vs sum ###########################
################################################################################################

dataSumAverageDoubletFinder = dataSumAverage %>% filter(condition == "DoubletFinder") %>% select(-X,-TNR,-auroc) %>% spread(type,auprc) %>% mutate(diff = sum - average)
dataSumAverageHybrid = dataSumAverage %>% filter(condition == "hybrid") %>% select(-X,-TNR,-auroc) %>% spread(type,auprc) %>% mutate(diff = sum - average)
dataSumAverageScDblFinder = dataSumAverage %>% filter(condition == "scDblFinder") %>% select(-X,-TNR,-auroc) %>% spread(type,auprc) %>% mutate(diff = sum - average)
dataSumAverageScrublet = dataSumAverage %>% filter(condition == "Scrublet") %>% select(-X,-TNR,-auroc) %>% spread(type,auprc) %>% mutate(diff = sum - average)

wilcox.test(dataSumAverageDoubletFinder$average, dataSumAverageDoubletFinder$sum, alternative = "two.sided", paired = TRUE, exact = TRUE)
wilcox.test(dataSumAverageHybrid$average, dataSumAverageHybrid$sum, paired = TRUE, alternative = "two.sided", exact = TRUE)
wilcox.test(dataSumAverageScDblFinder$average, dataSumAverageScDblFinder$sum, paired = TRUE, alternative = "two.sided", exact = TRUE)
wilcox.test(dataSumAverageScrublet$average, dataSumAverageScrublet$sum, paired = TRUE, alternative = "two.sided", exact = TRUE)


dataSumAverageDoubletFinder = dataSumAverage %>% filter(condition == "DoubletFinder") %>% select(-X,-TNR,-auprc) %>% spread(type,auroc) %>% mutate(diff = sum - average)
dataSumAverageHybrid = dataSumAverage %>% filter(condition == "hybrid") %>% select(-X,-TNR,-auprc) %>% spread(type,auroc) %>% mutate(diff = sum - average)
dataSumAverageScDblFinder = dataSumAverage %>% filter(condition == "scDblFinder") %>% select(-X,-TNR,-auprc) %>% spread(type,auroc) %>% mutate(diff = sum - average)
dataSumAverageScrublet = dataSumAverage %>% filter(condition == "Scrublet") %>% select(-X,-TNR,-auprc) %>% spread(type,auroc) %>% mutate(diff = sum - average)

wilcox.test(dataSumAverageDoubletFinder$average, dataSumAverageDoubletFinder$sum, alternative = "two.sided", paired = TRUE, exact = TRUE) #0.001677
wilcox.test(dataSumAverageHybrid$average, dataSumAverageHybrid$sum, paired = TRUE, alternative = "two.sided", exact = TRUE) # < 2.2e-16
wilcox.test(dataSumAverageScDblFinder$average, dataSumAverageScDblFinder$sum, paired = TRUE, alternative = "two.sided", exact = TRUE) # < 2.2e-16
wilcox.test(dataSumAverageScrublet$average, dataSumAverageScrublet$sum, paired = TRUE, alternative = "two.sided", exact = TRUE) # < 1.1e-12


dataSumAverageDoubletFinder = dataSumAverage %>% filter(condition == "DoubletFinder") %>% select(-X,-auroc,-auprc) %>% spread(type,TNR) %>% mutate(diff = sum - average)
dataSumAverageHybrid = dataSumAverage %>% filter(condition == "hybrid") %>% select(-X,-auroc,-auprc) %>% spread(type,TNR) %>% mutate(diff = sum - average)
dataSumAverageScDblFinder = dataSumAverage %>% filter(condition == "scDblFinder") %>% select(-X,-auroc,-auprc) %>% spread(type,TNR) %>% mutate(diff = sum - average)
dataSumAverageScrublet = dataSumAverage %>% filter(condition == "Scrublet") %>% select(-X,-auroc,-auprc) %>% spread(type,TNR) %>% mutate(diff = sum - average)

wilcox.test(dataSumAverageDoubletFinder$average, dataSumAverageDoubletFinder$sum, alternative = "two.sided", paired = TRUE, exact = TRUE) #0.003959
wilcox.test(dataSumAverageHybrid$average, dataSumAverageHybrid$sum, paired = TRUE, alternative = "two.sided", exact = TRUE) # < 2.2e-16
wilcox.test(dataSumAverageScDblFinder$average, dataSumAverageScDblFinder$sum, paired = TRUE, alternative = "two.sided", exact = TRUE) # < 2.2e-16
wilcox.test(dataSumAverageScrublet$average, dataSumAverageScrublet$sum, paired = TRUE, alternative = "two.sided", exact = TRUE) # 2.52e-05
