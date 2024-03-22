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

output_dir <- "/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp"
SPLINTR_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/SPLINTR/data/ten.rds")
LARRY_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/LARRY/data/ten.rds")
TREX_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/TREX/data/ten.rds")
FM03_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/FM03/data/ten.rds")
FM02_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/FM02/data/ten.rds")
FM01_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/FM01/data/ten.rds")
FM04_data <- readRDS("/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset/FM04/data/ten.rds")
FM05_data <- readRDS("/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset/FM05/data/ten.rds")
FM06_data <- readRDS("/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset/FM06/data/ten.rds")
FM08_data <- readRDS("/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset/FM08/data/ten.rds")
watermelon_data <- readRDS("/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset/watermelon/data/ten.rds")
non_cancer_data <- readRDS("/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset/non_cancer/data/ten.rds")
ClonMapper_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/ClonMapper/data/ten.rds")
cellTag_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/cellTag/data/ten.rds")
# Biorxiv_data <- readRDS("/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/Biorxiv/data/ten.rds")



data_list <- list(
  "LARRY"=LARRY_data,
  "SPLINTR"=SPLINTR_data,
  "FM01"=FM01_data,
  "FM04"=FM04_data,
  "FM05"=FM05_data,
  "FM06"=FM06_data,
  "FM08"=FM08_data,
  "TREX"=TREX_data,
  "FM03"=FM03_data,
  "FM02"=FM02_data,
  "watermelon"=watermelon_data,
  "non_cancer"=non_cancer_data,
  "ClonMapper"=ClonMapper_data,
  "cellTag"=cellTag_data
  # "Biorxiv"=Biorxiv_data
)

print("Loading finished!")

total_doublets <- 50
total_singlets <- 450
result_ls <- NULL

for(cur_dataset in names(data_list)){
  # if(cur_dataset == "SPLINTR" | cur_dataset == "LARRY"){
  #   total_doublets <- 300
  #   total_singlets <- 2700
  # }else{

  # }
  cur_s <- data_list[[cur_dataset]]
  cur_s <- JoinLayers(cur_s)
  cur_s@meta.data[["barcode"]] <- cur_s@meta.data|>rownames()
  cluster_count <- cur_s$RNA_snn_res.0.5|>table()|>as.data.frame()
  cluster_count <- cluster_count[order(cluster_count[,2], decreasing = T), ]

  adj_cluster <- cluster_count[1,1]
  adj_cell_id_doublet_pool <- cur_s@meta.data|>dplyr::filter(RNA_snn_res.0.5 == adj_cluster & label == "doublet")|>
    rownames()
  adj_cell_id_singlet_pool <- cur_s@meta.data|>dplyr::filter(RNA_snn_res.0.5 == adj_cluster & label == "singlet")|>
    rownames()

  three_cluster <- cluster_count[1:3,1]
  three_cluster_cell_id_doublet_pool <- cur_s@meta.data|>dplyr::filter(RNA_snn_res.0.5 %in% three_cluster & label == "doublet")|>
    rownames()
  three_cluster_cell_id_singlet_pool <- cur_s@meta.data|>dplyr::filter(RNA_snn_res.0.5 %in% three_cluster & label == "singlet")|>
    rownames()

  all_cluster_cell_id_doublet_pool <- cur_s@meta.data|>dplyr::filter(label == "doublet")|>rownames()
  all_cluster_cell_id_singlet_pool <- cur_s@meta.data|>dplyr::filter(label == "singlet")|>rownames()

  adj_test_1_cell_id <- c(sample(adj_cell_id_doublet_pool, total_doublets),
                          sample(adj_cell_id_singlet_pool, total_singlets))
  adj_test_2_cell_id <- c(sample(adj_cell_id_doublet_pool, total_doublets),
                          sample(adj_cell_id_singlet_pool, total_singlets))
  adj_test_3_cell_id <- c(sample(adj_cell_id_doublet_pool, total_doublets),
                          sample(adj_cell_id_singlet_pool, total_singlets))

  three_cluster_1_cell_id <- c(sample(three_cluster_cell_id_doublet_pool, total_doublets),
                               sample(three_cluster_cell_id_singlet_pool, total_singlets))
  three_cluster_2_cell_id <- c(sample(three_cluster_cell_id_doublet_pool, total_doublets),
                               sample(three_cluster_cell_id_singlet_pool, total_singlets))
  three_cluster_3_cell_id <- c(sample(three_cluster_cell_id_doublet_pool, total_doublets),
                               sample(three_cluster_cell_id_singlet_pool, total_singlets))

  all_cluster_1_cell_id <- c(sample(all_cluster_cell_id_doublet_pool, total_doublets),
                             sample(all_cluster_cell_id_singlet_pool, total_singlets))
  all_cluster_2_cell_id <- c(sample(all_cluster_cell_id_doublet_pool, total_doublets),
                             sample(all_cluster_cell_id_singlet_pool, total_singlets))
  all_cluster_3_cell_id <- c(sample(all_cluster_cell_id_doublet_pool, total_doublets),
                             sample(all_cluster_cell_id_singlet_pool, total_singlets))

  cur_dataset_tests <- list(
    "adj_test_1" = adj_test_1_cell_id,
    "adj_test_2" = adj_test_2_cell_id,
    "adj_test_3" = adj_test_3_cell_id,
    "three_cluster_1" = three_cluster_1_cell_id,
    "three_cluster_2" = three_cluster_2_cell_id,
    "three_cluster_3" = three_cluster_3_cell_id,
    "all_cluster_1" = all_cluster_1_cell_id,
    "all_cluster_2" = all_cluster_2_cell_id,
    "all_cluster_3" = all_cluster_3_cell_id
  )

  for(cur_test_name in names(cur_dataset_tests)){
    cell_ids <- cur_dataset_tests[[cur_test_name]]
    cur_obj <- subset(cur_s, subset = barcode %in% cell_ids)
    cur_obj_pc <- cur_obj@reductions$pca@cell.embeddings[,1:30]
    distances_cur_obj_pc_mean <- dist(cur_obj_pc)|>mean()
    print(cur_obj$label|>table())
    cur_obj_id <- glue("{cur_dataset}__{cur_test_name}")
    print(glue("INFO: {cur_obj_id}"))
    cur_type <- sub("(.*)_.+$", "\\1", cur_test_name)
    cur_obj <- JoinLayers(cur_obj)
    cur_sce <- as.SingleCellExperiment(cur_obj)
    sce <- scDblFinder::scDblFinder(cur_sce, dbr = 0.1)
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
      "type"=cur_type,
      "dataset_id"=cur_dataset,
      "mean_dist_pc"=distances_cur_obj_pc_mean
    )
    unique_id <- glue("doublet_cell__{cur_obj_id}")
    result_ls[[unique_id]] <- output
    obj_out <- glue("/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp/objects/{unique_id}.rds")
    saveRDS(sce, obj_out)

    doublet.seurat <- NormalizeData(cur_obj)
    doublet.seurat <- ScaleData(doublet.seurat)
    doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
    doublet.seurat <- RunPCA(doublet.seurat)

    sweep.res.doublet <- paramSweep(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 8)
    sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
    bcmvn.doublet <- find.pK(sweep.stats.doublet)
    pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
    # nExp_poi <- sum(doublet.seurat[["label"]] == "doublet")
    nExp_poi<-round(dim(doublet.seurat)[2]*0.1)
    doublet.seurat <- doubletFinder(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)

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
      "type"=cur_type,
      "dataset_id"=cur_dataset,
      "mean_dist_pc"=distances_cur_obj_pc_mean
    )
    unique_id <- glue("doublet_finder__{cur_obj_id}")
    result_ls[[unique_id]] <- output

    obj_out <- glue("/projects/p31666/zzhang/doublet-bchmk/final/cell_location_comp/objects/{unique_id}.rds")
    saveRDS(doublet.seurat, obj_out)
  }
}

result_df <- result_ls|>rbindlist()|>as.data.frame()
out_file <- glue("{output_dir}/all.tsv")
write.table(result_df, file = out_file, quote = FALSE, sep="\t",row.names = FALSE)