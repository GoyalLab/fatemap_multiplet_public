run_scDblFinder_on_dataset <- function(data_dir, METHOD_ID, dbl_exp, dbl_act, num_to_test=5,
                                        save_root = "/projects/b1042/GoyalLab/zzhang/doublet_objects/scDblFinder") {
  print(glue("Working on dataset: {data_dir}"))
  dir_10X <- glue("{data_dir}/10X")
  dataset_id <- strsplit(data_dir, "/")[[1]][[length(strsplit(data_dir, "/")[[1]])]]
  sample_dirs <- list.dirs(path = dir_10X, full.names = TRUE, recursive = F)
  stats <- list()
  out_df <- NULL
  for (cur_sample_dir in sample_dirs) {
    cur_sample_id <- strsplit(cur_sample_dir, "/")[[1]][[length(strsplit(cur_sample_dir, "/")[[1]])]]
    cur_data_sample_id <- glue("{dataset_id}_{cur_sample_id}")
    cur_save_dir <- glue("{save_root}/act__{dbl_act}")
    if (!file.exists(cur_save_dir)) {
      dir.create(cur_save_dir, recursive = TRUE)
    }
    cur_df <- run_bchmk(cur_sample_dir = cur_sample_dir, dbl_act = dbl_act,
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
    cur_file <- glue("/projects/b1042/GoyalLab/zzhang/doublet_objects_benchmarking_data/scDblFinder/{cur_data_id}.rds")
    cur_res <- run_scDblFinder(sce, cur_exp_dbl, doublet_obj_file = cur_file)
    df_ls <- rbindlist(list(df_ls, cur_res))
  }
  all_stats <- data.frame(df_ls)
  all_stats <- cbind(ID = rep(cur_data_id, nrow(all_stats)), all_stats)
  all_stats
}


run_bchmk <- function(cur_sample_dir, dbl_act, data_sample_ID, num_to_test, dbl_exp, save_dir) {
  df_ls <- NULL
  for (idx in 1:num_to_test) {
    sce <- load_fatemap_into_sce(cur_sample_dir, doublets_pct=dbl_act)
    prefix1 <- basename(cur_sample_dir)
    prefix2 <- tail(strsplit(cur_sample_dir, "/")[[1]], 3) |>
      head(1)
    doublet_object_file <- glue("{save_dir}/{prefix2}___{prefix1}___exp_{dbl_exp}__act_{dbl_act}.rds")
    cur_res <- run_scDblFinder(sce, dbl_exp, doublet_object_file)
    df_ls <- rbindlist(list(df_ls, cur_res))

  }
  all_stats <- data.frame(df_ls)
  all_stats <- cbind(ID = rep(data_sample_ID, nrow(all_stats)), all_stats)
  all_stats
}

run_scDblFinder <- function(sce, dbl_exp, doublet_obj_file) {

  # print(glue("INFO: Current dataset doublet rate is {cur_doublet_rate}"))
  sce <- scDblFinder::scDblFinder(sce, dbr = dbl_exp)

  # save
  if(is.null(doublet_obj_file) == FALSE){
    saveRDS(sce, doublet_obj_file)
  }

  score <- sce[["scDblFinder.score"]]
  fg <- score[sce$label=="doublet"]
  bg <- score[sce$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc
  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr
  )
}

load_fatemap_into_sce <- function(data_dir, doublets_pct = 0.08) {
  doublet.seurat <- load_fatemap_into_seurat(data_dir = data_dir, doublets_pct = doublets_pct)
  sce <- as.SingleCellExperiment(doublet.seurat)
  sce
}

load_fatemap_into_seurat <- function(data_dir, doublets_pct = 0.08, dbl_method = "sum") {
  expression_matrix <- Read10X(data.dir = data_dir)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  singlets_file <- glue("{data_dir}/singlets_all.txt")

  # key file for simulaitng doublets
  singlets_pair_for_doublet_simulation_file <- list.files(data_dir, pattern = "corrected_singlet_pairs.csv")

  # extract the singlets
  singlets_pair_for_doublet_simulation_file<-glue("{data_dir}/{singlets_pair_for_doublet_simulation_file}")

  # check if "-1" is a part of 10X barcodes
  if(!grepl("watermelon", data_dir)){
    if(grepl("-1", seurat_object|>colnames()|>head(1))){
      singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
        as.vector()|>
        paste("-1", sep = "")
    }else{
      singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
        as.vector()|>
        str_remove("-1")
    }
  }else{
    plus_one_str <- sub(".*(-\\d+)$", "\\1", seurat_object|>colnames()|>head(1))
    singlets <- read.delim(singlets_file, header = F, sep = "\n")[, 1]|>
      as.vector()|>
      paste(plus_one_str, sep = "")
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
  if(dbl_method == "avg"){
    print("INFO: avg doublets!")
    avg_assay_ori <- (assay_for_simulated_doublets[, real_cells_pt1] + assay_for_simulated_doublets[, real_cells_pt2]) / 2
  }else{
    print("INFO: summation doublets!")
    avg_assay_ori <- assay_for_simulated_doublets[, real_cells_pt1] + assay_for_simulated_doublets[, real_cells_pt2]
  }

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
