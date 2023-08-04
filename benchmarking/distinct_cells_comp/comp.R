library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(PRROC)
library(pbapply)
library(scDblFinder)
require(glue)
require(purrr)
require(data.table)
library(stringr)
require(parallel)
library(DoubletFinder)
library(ggplot2)
library(tidyverse)


options(error=traceback)

theme_Publication_blank <- function(base_size=12, base_family="", lgd_position="bottom") { #12 For ALDR paper
  require(grid)
  require(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.border = element_rect(colour = NA, fill = "transparent"),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,margin=margin(0,10,0,0)),
            axis.title.x = element_text(margin=margin(10,0,0,0)),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(size = 0.3),
            axis.line.x = element_line(size = 0.3, linetype = "solid", colour = "black"),
            axis.line.y = element_line(size = 0.3, linetype = "solid", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA, fill="transparent"),
            legend.position = lgd_position,
            #legend.direction = "horizontal",
            #legend.box = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(10, "pt"),
            #legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#d8d8d8",fill="#d8d8d8")
            # strip.text = element_text(face="bold")
  ))

}

FM04_s <- readRDS("/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp/FM04_integrated.rds")
DefaultAssay(FM05_s)<-"RNA"
FM03_s <- readRDS("/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp/FM03_integrated.rds")
DefaultAssay(FM01_s)<-"RNA"
FM02_s <- readRDS("/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp/FM02_integrated.rds")
DefaultAssay(FM02_s)<-"RNA"
output_dir <- "/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp"

p1 <- DimPlot(FM02_s, reduction = "umap", repel = TRUE, label = T, group.by = "integrated_snn_res.0.5") + ggtitle("FM02")
p2 <- DimPlot(FM03_s, reduction = "umap", repel = TRUE, label = T, group.by = "integrated_snn_res.0.5") + ggtitle("FM03")
p3 <- DimPlot(FM04_s, reduction = "umap", repel = TRUE, label = T, group.by = "integrated_snn_res.0.8") + ggtitle("FM04")

out_file <- glue("{output_dir}/umaps_combined.svg")
svg(out_file, width = 15, height = 6)
plot(p1 + p2 + p3)
dev.off()
#distant
FM03_clus_4_5 <- subset(FM03_s, integrated_snn_res.0.5==4 | integrated_snn_res.0.5 ==5)
FM03_clus_3_6 <- subset(FM03_s, integrated_snn_res.0.5==3 | integrated_snn_res.0.5 ==6)
FM02_clus_0_4 <- subset(FM02_s, integrated_snn_res.0.5==0 | integrated_snn_res.0.5 ==4)
FM02_clus_6_7 <- subset(FM02_s, integrated_snn_res.0.5==6 | integrated_snn_res.0.5 ==7)
FM04_clus_1_3 <- subset(FM04_s, integrated_snn_res.0.8==1 | integrated_snn_res.0.8 ==3)
FM04_clus_8_4 <- subset(FM04_s, integrated_snn_res.0.8==8 | integrated_snn_res.0.8 ==4)

distant_tests <- c(
  "FM03_clus_4_5" = FM03_clus_4_5,
  "FM03_clus_3_6" = FM03_clus_3_6,
  "FM02_clus_0_4" = FM02_clus_0_4,
  "FM02_clus_6_7" = FM02_clus_6_7,
  "FM04_clus_1_3" = FM04_clus_1_3,
  "FM04_clus_8_4" = FM04_clus_8_4
)

#adjacent
FM03_clus_3_4 <- subset(FM03_s, integrated_snn_res.0.5==3 | integrated_snn_res.0.5 ==4)
FM03_clus_5_6 <- subset(FM03_s, integrated_snn_res.0.5==5 | integrated_snn_res.0.5 ==6)
FM02_clus_4_6 <- subset(FM02_s, integrated_snn_res.0.5==4 | integrated_snn_res.0.5 ==6)
FM02_clus_0_2 <- subset(FM02_s, integrated_snn_res.0.5==0 | integrated_snn_res.0.5 ==2)
FM04_clus_1_5 <- subset(FM04_s, integrated_snn_res.0.8==1 | integrated_snn_res.0.8 ==5)
FM04_clus_3_6 <- subset(FM04_s, integrated_snn_res.0.8==3 | integrated_snn_res.0.8 ==6)

adjacent_tests <- c(
  "FM03_clus_3_4" = FM03_clus_3_4,
  "FM03_clus_5_6" = FM03_clus_5_6,
  "FM02_clus_4_6" = FM02_clus_4_6,
  "FM02_clus_0_2" = FM02_clus_0_2,
  "FM04_clus_1_5" = FM04_clus_1_5,
  "FM04_clus_3_6" = FM04_clus_3_6
)

# baseline
set.seed(2022)
FM03_rand_1 <- FM03_s[, sample(colnames(FM03_s), size = 2000, replace=F)]
set.seed(2023)
FM02_rand_1 <- FM02_s[, sample(colnames(FM02_s), size = 2000, replace=F)]
set.seed(2024)
FM04_rand_1 <- FM04_s[, sample(colnames(FM04_s), size = 2000, replace=F)]
set.seed(2025)
FM03_rand_2 <- FM03_s[, sample(colnames(FM03_s), size = 2000, replace=F)]
set.seed(2026)
FM02_rand_2 <- FM02_s[, sample(colnames(FM02_s), size = 2000, replace=F)]
set.seed(2027)
FM04_rand_2 <- FM04_s[, sample(colnames(FM04_s), size = 2000, replace=F)]
baseline_tests <- c(
  "FM03_rand_1" = FM03_rand_1,
  "FM03_rand_2" = FM03_rand_2,
  "FM02_rand_1" = FM02_rand_1,
  "FM02_rand_2" = FM02_rand_2,
  "FM04_rand_1" = FM04_rand_1,
  "FM04_rand_2" = FM04_rand_2
)



saveRDS(object = baseline_tests, file = glue("{output_dir}/baseline_tests.rds"))
saveRDS(object = distant_tests, file = glue("{output_dir}/distant_tests.rds"))
saveRDS(object = adjacent_tests, file = glue("{output_dir}/adjacent_tests.rds"))

result_ls <- NULL

for(cur_obj_id in names(baseline_tests)){
  print(glue("INFO: Currently working on {cur_obj_id}"))
  cur_dataset_id <- str_split(cur_obj_id, "_")[[1]][1]
  cur_obj <- baseline_tests[[cur_obj_id]]
  cur_sce <- as.SingleCellExperiment(cur_obj)
  sce <- scDblFinder::scDblFinder(cur_sce, dbr = 0.08)
  score <- sce[["scDblFinder.score"]]
  fg <- score[sce$label=="doublet"]
  bg <- score[sce$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc
  print(glue("ROC: {cur_roc}! PR: {cur_pr}"))
  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr,
    "method"="doublet_cell",
    "type"="baseline",
    "dataset_id"=cur_dataset_id
  )|>
    data.frame()
  result_ls <- rbindlist(list(result_ls, output))

  doublet.seurat <- NormalizeData(cur_obj)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)

  sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 8)
  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
  # nExp_poi <- sum(doublet.seurat[["label"]] == "doublet")
  nExp_poi<-round(dim(doublet.seurat)[2]*0.08)
  doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)

  attribute <- paste('pANN', 0.25, pK, nExp_poi, sep = '_')
  classification_label<-paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')

  score <- doublet.seurat@meta.data[[attribute]]
  fg <- score[doublet.seurat@meta.data$label=="doublet"]
  bg <- score[doublet.seurat@meta.data$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc

  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr,
    "method"="doublet_finder",
    "type"="baseline",
    "dataset_id"=cur_dataset_id
  )|>
    data.frame()
  result_ls <- rbindlist(list(result_ls, output))
}

for(cur_obj_id in names(adjacent_tests)){
  print(glue("INFO: Currently working on {cur_obj_id}"))
  cur_dataset_id <- str_split(cur_obj_id, "_")[[1]][1]
  cur_obj <- adjacent_tests[[cur_obj_id]]
  cur_sce <- as.SingleCellExperiment(cur_obj)
  sce <- scDblFinder::scDblFinder(cur_sce, dbr = 0.08)
  score <- sce[["scDblFinder.score"]]
  fg <- score[sce$label=="doublet"]
  bg <- score[sce$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc
  print(glue("ROC: {cur_roc}! PR: {cur_pr}"))
  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr,
    "method"="doublet_cell",
    "type"="adjacent",
    "dataset_id"=cur_dataset_id
  )|>
    data.frame()
  result_ls <- rbindlist(list(result_ls, output))

  doublet.seurat <- NormalizeData(cur_obj)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)

  sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 8)
  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
  # nExp_poi <- sum(doublet.seurat[["label"]] == "doublet")
  nExp_poi<-round(dim(doublet.seurat)[2]*0.08)
  doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)

  attribute <- paste('pANN', 0.25, pK, nExp_poi, sep = '_')
  classification_label<-paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')

  score <- doublet.seurat@meta.data[[attribute]]
  fg <- score[doublet.seurat@meta.data$label=="doublet"]
  bg <- score[doublet.seurat@meta.data$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc

  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr,
    "method"="doublet_finder",
    "type"="adjacent",
    "dataset_id"=cur_dataset_id
  )|>
    data.frame()
  result_ls <- rbindlist(list(result_ls, output))
}


for(cur_obj_id in names(distant_tests)){
  print(glue("INFO: Currently working on {cur_obj_id}"))
  cur_dataset_id <- str_split(cur_obj_id, "_")[[1]][1]
  cur_obj <- distant_tests[[cur_obj_id]]
  cur_sce <- as.SingleCellExperiment(cur_obj)
  sce <- scDblFinder::scDblFinder(cur_sce, dbr = 0.08)
  score <- sce[["scDblFinder.score"]]
  fg <- score[sce$label=="doublet"]
  bg <- score[sce$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc
  print(glue("ROC: {cur_roc}! PR: {cur_pr}"))
  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr,
    "method"="doublet_cell",
    "type"="distant",
    "dataset_id"=cur_dataset_id
  )|>
    data.frame()
  result_ls <- rbindlist(list(result_ls, output))

  doublet.seurat <- NormalizeData(cur_obj)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)

  sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 8)
  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
  # nExp_poi <- sum(doublet.seurat[["label"]] == "doublet")
  nExp_poi<-round(dim(doublet.seurat)[2]*0.08)
  doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)

  attribute <- paste('pANN', 0.25, pK, nExp_poi, sep = '_')
  classification_label<-paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')

  score <- doublet.seurat@meta.data[[attribute]]
  fg <- score[doublet.seurat@meta.data$label=="doublet"]
  bg <- score[doublet.seurat@meta.data$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc

  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr,
    "method"="doublet_finder",
    "type"="distant",
    "dataset_id"=cur_dataset_id
  )|>
    data.frame()
  result_ls <- rbindlist(list(result_ls, output))
}

result_df <- result_ls|>as.data.frame()
out_file <- glue("{output_dir}/all.tsv")
write.table(result_df, file = out_file, quote = FALSE, sep="\t",row.names = FALSE)

# Generate plots
# Calculate mean ROC for each type, method and dataset
x_means_dataset <- result_df %>%
  group_by(type, method, dataset_id) %>%
  summarise(mean_roc = mean(roc), .groups = 'drop')

# Calculate mean ROC for each type and method
x_means_method <- x_means_dataset %>%
  group_by(type, method) %>%
  summarise(mean_roc = mean(mean_roc), .groups = 'drop')

# Create the scatter plot
p <- ggplot(result_df, aes(x = type, y = roc, color = method)) +
  geom_point(aes(shape = dataset_id), width = 0.2, size = 3) +  # Use geom_jitter to plot all data points and shape by dataset_id
  geom_line(data = x_means_method, aes(y = mean_roc, group = method)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Type", y = "ROC", color = "Method", shape = "Dataset ID")+
  ggtitle("ROC")+
  theme_Publication_blank(lgd_position = "right")

# 6 lines

x_means <- result_df %>%
  group_by(type, method, dataset_id) %>%
  summarise(mean_roc = mean(roc), .groups = 'drop')

# Create the scatter plot
p2 <- ggplot(result_df, aes(x = type, y = roc, color = method, shape = dataset_id)) +
  geom_point(width = 0.2, size = 3) +  # Use geom_jitter to plot all data points
  geom_line(data = x_means, aes(y = mean_roc, group = interaction(dataset_id, method))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Type", y = "ROC", color = "Method", shape = "Dataset ID")+
  ggtitle("ROC")+
  theme_Publication_blank(lgd_position = "right")
# Print the plot

out_file <- glue("{output_dir}/roc.svg")
svg(out_file, height=6, width=8)
print(p)
dev.off()


######            PR              ################

# Calculate mean pr for each type, method and dataset
x_means_dataset <- result_df %>%
  group_by(type, method, dataset_id) %>%
  summarise(mean_pr = mean(pr), .groups = 'drop')

# Calculate mean pr for each type and method
x_means_method <- x_means_dataset %>%
  group_by(type, method) %>%
  summarise(mean_pr = mean(mean_pr), .groups = 'drop')

# Create the scatter plot
p <- ggplot(result_df, aes(x = type, y = pr, color = method)) +
  geom_point(aes(shape = dataset_id), width = 0.2, size = 3) +  # Use geom_jitter to plot all data points and shape by dataset_id
  geom_line(data = x_means_method, aes(y = mean_pr, group = method)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Type", y = "pr", color = "Method", shape = "Dataset ID")+
  ggtitle("pr")+
  theme_Publication_blank(lgd_position = "right")

# 6 lines

x_means <- result_df %>%
  group_by(type, method, dataset_id) %>%
  summarise(mean_pr = mean(pr), .groups = 'drop')

# Create the scatter plot
p2 <- ggplot(result_df, aes(x = type, y = pr, color = method, shape = dataset_id)) +
  geom_point(width = 0.2, size = 3) +  # Use geom_jitter to plot all data points
  geom_line(data = x_means, aes(y = mean_pr, group = interaction(dataset_id, method))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Type", y = "pr", color = "Method", shape = "Dataset ID")+
  ggtitle("pr")+
  theme_Publication_blank(lgd_position = "right")
# Print the plot
out_file <- glue("{output_dir}/pr.svg")
svg(out_file, height=6, width=8)
print(p)
dev.off()


# UMAP for each test dataset

FM03_data <- c(
  "FM03_clus_4_5" = FM03_clus_4_5,
  "FM03_clus_3_6" = FM03_clus_3_6,
  "FM03_clus_3_4" = FM03_clus_3_4,
  "FM03_clus_5_6" = FM03_clus_5_6,
  "FM03_rand_1" = FM03_rand_1,
  "FM03_rand_2" = FM03_rand_2
)

for(cur_FM03_data_name in names(FM03_data)){
  cur_FM03_data <- FM03_data[[cur_FM03_data_name]]
  p <- DimPlot(FM03_s, cells.highlight = colnames(cur_FM03_data))+
    scale_color_manual(labels = c("Unselected", "Selected"), values = c("gray", "hotpink3"))+
    theme(legend.position = "none")+
    theme(axis.title = element_blank())+
    theme(axis.ticks = element_blank(), axis.text = element_blank())+
    theme(axis.line = element_blank())
  out_file <- glue("{output_dir}/{cur_FM03_data_name}.png")
  ggsave(filename = out_file, plot = p, width = 6, height = 6, dpi = 300)
}

FM02_data <- c(
  "FM02_clus_0_4" = FM02_clus_0_4,
  "FM02_clus_6_7" = FM02_clus_6_7,
  "FM02_clus_4_6" = FM02_clus_4_6,
  "FM02_clus_0_2" = FM02_clus_0_2,
  "FM02_rand_1" = FM02_rand_1,
  "FM02_rand_2" = FM02_rand_2
)

for(cur_FM02_data_name in names(FM02_data)){
  cur_FM02_data <- FM02_data[[cur_FM02_data_name]]
  p <- DimPlot(FM02_s, cells.highlight = colnames(cur_FM02_data))+
    scale_color_manual(labels = c("Unselected", "Selected"), values = c("gray", "hotpink3"))+
    theme(legend.position = "none")+
    theme(axis.title = element_blank())+
    theme(axis.ticks = element_blank(), axis.text = element_blank())+
    theme(axis.line = element_blank())
  out_file <- glue("{output_dir}/{cur_FM02_data_name}.png")
  ggsave(filename = out_file, plot = p, width = 6, height = 6, dpi = 300)
}

FM04_data <- c(
  "FM04_clus_1_3" = FM04_clus_1_3,
  "FM04_clus_8_4" = FM04_clus_8_4,
  "FM04_clus_1_5" = FM04_clus_1_5,
  "FM04_clus_3_6" = FM04_clus_3_6,
  "FM04_rand_1" = FM04_rand_1,
  "FM04_rand_2" = FM04_rand_2
)

for(cur_FM04_data_name in names(FM04_data)){
  cur_FM04_data <- FM04_data[[cur_FM04_data_name]]
  p <- DimPlot(FM04_s, cells.highlight = colnames(cur_FM04_data))+
    scale_color_manual(labels = c("Unselected", "Selected"), values = c("gray", "hotpink3"))+
    theme(legend.position = "none")+
    theme(axis.title = element_blank())+
    theme(axis.ticks = element_blank(), axis.text = element_blank())+
    theme(axis.line = element_blank())
  out_file <- glue("{output_dir}/{cur_FM04_data_name}.png")
  ggsave(filename = out_file, plot = p, width = 6, height = 6, dpi = 300)
}