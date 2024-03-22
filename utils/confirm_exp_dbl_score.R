library(Matrix)
library(Seurat)
library(scds)
library(data.table)
library(SingleCellExperiment)
library(PRROC)
require(DoubletFinder)
source("/projects/p31666/zzhang/doublet-bchmk/repo/benchmarking/scds/scds_utils.R")

root_dir <-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
METHOD_ID<-"hybrid"
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final"
root_bchmk_dataset<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
num_test <- 1
doublet_rates_to_test<-c(0.05, 0.08, 0.1, 0.15, 0.2, 0.25)
ACT_DBL <- 0.15
################################
#                              #
#             Main             #
#                              #
################################




run_hybrid_on_dataset <- function(data_dir, METHOD_ID, dbl_exp, dbl_act, num_to_test=5,
                                  save_root = "/projects/b1042/GoyalLab/zzhang/doublet_objects/hybrid") {
  print(glue("Working on dataset: {data_dir}"))
  dir_10X <- glue("{data_dir}/10X")
  dataset_id <- strsplit(data_dir, "/")[[1]][[length(strsplit(data_dir, "/")[[1]])]]
  sample_dirs <- list.dirs(path = dir_10X, full.names = TRUE, recursive = F)
  stats <- list()
  out_df <- NULL
  labels_prefix <- glue("{data_dir}/fatemapID/{dataset_id}")
  for (cur_sample_dir in sample_dirs) {
    cur_sample_id <- strsplit(cur_sample_dir, "/")[[1]][[length(strsplit(cur_sample_dir, "/")[[1]])]]
    cur_data_sample_id <- glue("{dataset_id}_{cur_sample_id}")
    cur_save_dir <- glue("{save_root}/act__{dbl_act}")
    if (!file.exists(cur_save_dir)) {
      dir.create(cur_save_dir, recursive = TRUE)
    }
    cur_df <- run_bchmk(cur_sample_dir = cur_sample_dir, labels_prefix = labels_prefix, dbl_act = dbl_act,
                        data_sample_ID = cur_data_sample_id, num_to_test = num_to_test, dbl_exp = dbl_exp,
                        save_dir = cur_save_dir)
    out_df <- rbindlist(list(out_df, cur_df))|>
      data.frame()
  }
  # out_df <- out_df[order(out_df$Expected_Doublet_Rate),]
  out_df
}

run_bchmk_original_paper <- function(cur_bchmk_dataset_mtx, num_to_test){
  df_ls <- NULL
  for(idx in 1:num_test){
    cur_mtx <- readMM(cur_bchmk_dataset_mtx)
    cur_data_id<-strsplit(cur_bchmk_dataset_mtx, "/")[[1]][[8]]
    cur_mtx_cell_id_file <- str_replace(cur_bchmk_dataset_mtx, "0.08.mtx", "cell_id_0.08.csv")
    cur_mtx_gene_id_file <- str_replace(cur_bchmk_dataset_mtx, "0.08.mtx", "gene_id_0.08.csv")
    cur_mtx_label <- str_replace(cur_bchmk_dataset_mtx, "0.08.mtx", "labels_0.08.csv")
    cur_mtx_cell_id <- read.table(cur_mtx_cell_id_file)[["V1"]]
    cur_mtx_gene_id <- read.table(cur_mtx_gene_id_file)[["V1"]]
    cur_mtx_label <- read.table(cur_mtx_label)[["V1"]]
    rownames(cur_mtx) <- cur_mtx_gene_id
    colnames(cur_mtx) <- cur_mtx_cell_id
    cur_mtx_label <- as.data.frame(cur_mtx_label)
    rownames(cur_mtx_label) <- cur_mtx_cell_id
    
    doublet.seurat <- CreateSeuratObject(counts = cur_mtx, min.cells = 1)
    doublet.seurat <- AddMetaData(doublet.seurat, col.name="label", cur_mtx_label)
    cur_exp_dbl<-sum(doublet.seurat[["label"]] == "doublet")/ncol(doublet.seurat)
    print(glue("INFO: Current benchmark dataset doublet rate is: {cur_exp_dbl}"))
    sce <- as.SingleCellExperiment(doublet.seurat)
    cur_res <- run_cxds(sce, exp_dbl = cur_exp_dbl, doublet_obj_file = NULL)
    df_ls <- rbindlist(list(df_ls, cur_res))
  }
  all_stats <- data.frame(df_ls)
  all_stats <- cbind(ID = rep(cur_data_id, nrow(all_stats)), all_stats)
  all_stats
}

run_bchmk <- function(cur_sample_dir, labels_prefix, dbl_act, data_sample_ID, num_to_test, dbl_exp, save_dir) {
  df_ls <- NULL
  for (idx in 1:num_to_test) {
    sce <- load_fatemap_into_sce(cur_sample_dir, labels_prefix, doublets_pct=dbl_act)
    prefix1 <- basename(cur_sample_dir)
    prefix2 <- basename(labels_prefix)
    doublet_object_file <- glue("{save_dir}/{prefix2}___{prefix1}___exp_{dbl_exp}__act_{dbl_act}.rds")
    cur_res <- run_cxds(sce, exp_dbl=dbl_exp, doublet_obj_file=doublet_object_file)
    df_ls <- rbindlist(list(df_ls, cur_res))
    
  }
  all_stats <- data.frame(df_ls)
  all_stats <- cbind(ID = rep(data_sample_ID, nrow(all_stats)), all_stats)
  all_stats
}

run_cxds <- function(sce, exp_dbl, doublet_obj_file) {
  n_exp_dbl <- round(dim(sce)[2]*exp_dbl)
  print(glue("Expected doublet rate: {exp_dbl}. Expected doublets: {n_exp_dbl}"))
  sce <- cxds_bcds_hybrid(sce, estNdbl = n_exp_dbl, force=T)
  
  # save
  if(is.null(doublet_obj_file) == FALSE){
    saveRDS(sce, doublet_obj_file)
  }
  
  score <- colData(sce)
  fg <- score[score$label=="doublet",][["hybrid_score"]]
  bg <- score[score$label=="singlet",][["hybrid_score"]]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc
  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr
  )
}

load_fatemap_into_sce <- function(data_dir, cell_labels_prefix, doublets_pct = 0.08) {
  doublet.seurat <- load_fatemap_into_seurat(data_dir = data_dir, cell_labels_prefix =
                                               cell_labels_prefix, doublets_pct = doublets_pct)
  sce <- as.SingleCellExperiment(doublet.seurat)
  sce
}

load_fatemap_into_seurat <- function(data_dir, cell_labels_prefix, doublets_pct = 0.08) {
  expression_matrix <- Read10X(data.dir = data_dir)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  singlets_file <- glue("{cell_labels_prefix}_singlets_all.txt")
  
  # key file for simulaitng doublets
  singlets_pair_for_doublet_simulation_file <- list.files(data_dir, pattern = "corrected_singlet_pairs.csv")
  
  # extract the singlets
  singlets_pair_for_doublet_simulation_file<-glue("{data_dir}/{singlets_pair_for_doublet_simulation_file}")
  # check if "-1" is a part of 10X barcodes
  if(grepl("-1", seurat_object|>colnames()|>head(1))){
    singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
      as.vector()|>
      paste("-1", sep = "")
  }else{
    singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
      as.vector()|>
      str_remove("-1")
  }
  
  # extract the singlet pairs
  singlet_pair_df<-read.csv(singlets_pair_for_doublet_simulation_file, header = F)|>
    #dplyr::mutate_each(funs = function(x) str_remove_all(x," "))
    dplyr::mutate(across(everything(), ~str_remove_all(.x, " ")))
  # remove combinations that did not pass QC in seurat
  all_seurat_cells<-seurat_object@meta.data|>row.names()
  
  # take average of expression to produce doublets
  real_cells_pt1 <- singlet_pair_df$V1
  real_cells_pt2 <- singlet_pair_df$V2
  assay_for_simulated_doublets <- seurat_object[["RNA"]]$counts[, c(real_cells_pt1, real_cells_pt2)]
  avg_assay_ori <- (assay_for_simulated_doublets[, real_cells_pt1] + assay_for_simulated_doublets[, real_cells_pt2]) / 2
  
  # create artificial doublet ID
  singlets_in_seurat <- singlets[singlets %in% colnames(seurat_object)]
  doublet_id<-paste(singlet_pair_df$V1,singlet_pair_df$V2, sep = "--")
  colnames(avg_assay_ori)<-doublet_id
  num_doublets<-round(length(singlets_in_seurat)*doublets_pct/(1+doublets_pct)) #derived from the case where 2 unique singlets make up 1 doublet, num_doublets = doublets_pct*(num_singlets-num_doublets)
  avg_assay <- avg_assay_ori[,1:num_doublets]
  doublet_object<-CreateSeuratObject(counts = avg_assay)
  
  
  singlets_used_for_doublets <- unique(unlist(strsplit(colnames(avg_assay), "--")))
  singlets_to_keep <- setdiff(colnames(seurat_object), singlets_used_for_doublets)
  singlets_to_keep_in_seurat <- intersect(singlets_in_seurat, singlets_to_keep)
  if (length(singlets_to_keep_in_seurat) > (num_doublets*(1-doublets_pct)/doublets_pct)) {
    singlets_trimmed <- sample(singlets_to_keep_in_seurat, num_doublets*(1-doublets_pct)/doublets_pct)
  } else {
    num_doublets <- round(length(singlets_in_seurat)*doublets_pct/(1-doublets_pct)) #calculate the number of doublets, no assumptions about removal being made
    avg_assay_new<-avg_assay_ori[,1:num_doublets]
    singlets_used_for_doublets <- unique(unlist(strsplit(colnames(avg_assay_new), "--")))
    singlets_to_keep <- setdiff(colnames(seurat_object), singlets_used_for_doublets)
    singlets_to_keep_in_seurat <- intersect(singlets_in_seurat, singlets_to_keep) # remove singlets used to make doublets
    num_doublets_adjusted <- round((length(singlets_to_keep_in_seurat)*doublets_pct) / (1-doublets_pct)) # adjust the number of doublets based on the number of singlets removed
    avg_assay_adjusted <- avg_assay_new[, 1:num_doublets_adjusted] # subset the avg_assay_new to be as long as the adjusted doublet count
    doublet_object <- CreateSeuratObject(counts = avg_assay_adjusted) #create doublet object
    singlets_trimmed <- singlets_to_keep_in_seurat #no more singlets removed after the ones used to make doublets are removed
  }
  singlet_object <- seurat_object[, singlets_trimmed]
  
  # tag doublet labels back to the doublet object
  label<-rep("doublet",ncol(doublet_object))|>
    as.data.frame()
  rownames(label)<-colnames(doublet_object)
  doublet_object<-AddMetaData(doublet_object, label, col.name = "label")
  
  # merge the two objects
  label<-rep("singlet",ncol(singlet_object))|>
    as.data.frame()
  rownames(label)<-colnames(singlet_object)
  singlet_object<-AddMetaData(singlet_object, label, col.name = "label")
  doublet.seurat<-merge(singlet_object, y = doublet_object, project = "combined_doublet_singlet")
  doublet_count <- sum(doublet.seurat@meta.data$label == "doublet")
  print(doublet_count)
  singlet_count <- sum(doublet.seurat@meta.data$label == "singlet")
  print(singlet_count)
  # Seurat V5 modifications
  t <- JoinLayers(doublet.seurat)
  doublet.seurat[["RNA"]] <- t[["RNA"]]
  doublet.seurat
}











main<-function (root_dir, METHOD_ID, out_dir){
  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL
  
  
  for(expected_doublet_rate in doublet_rates_to_test){
      cur_dataset_dir <- "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM01"
      print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{ACT_DBL}"))
      cur_out_df<-run_hybrid_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
                                        dbl_exp = expected_doublet_rate, dbl_act = ACT_DBL)
      cur_out_df$dbl_exp <- expected_doublet_rate
      cur_out_df$dbl_act <- ACT_DBL
      out_df<-rbindlist(list(out_df, cur_out_df))|>
        data.frame()
      
  }
  
  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_{ACT_DBL}.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}


main(root_dir, METHOD_ID, out_dir)