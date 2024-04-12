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
  singlet_count <- sum(doublet.seurat@meta.data$label == "singlet")
  # Seurat V5 modifications
  t <- JoinLayers(doublet.seurat)
  doublet.seurat[["RNA"]] <- t[["RNA"]]
  doublet.seurat
}

run_doublet_finder<-function(doublet.seurat, exp_dbl, doublet_obj_file){

  doublet.seurat <- NormalizeData(doublet.seurat)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)

  sweep.res.doublet <- paramSweep(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 8)
  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
  # nExp_poi <- sum(doublet.seurat[["label"]] == "doublet")
  nExp_poi<-round(dim(doublet.seurat)[2]*exp_dbl)
  doublet.seurat <- doubletFinder(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)

  attribute <- paste('pANN', 0.25, pK, nExp_poi, sep = '_')
  classification_label<-paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')

  # save
  if(is.null(doublet_obj_file) == FALSE){
    saveRDS(doublet.seurat, doublet_obj_file)
  }

  score <- doublet.seurat@meta.data[[attribute]]
  fg <- score[doublet.seurat@meta.data$label=="doublet"]
  bg <- score[doublet.seurat@meta.data$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc
  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr
  )
  output
}
run_bchmk<-function(data_dir, data_sample_ID, num_to_test, exp_dbl, act_dbl, save_dir){
  df_ls<-NULL
  for(idx in 1:num_to_test){
    seurat_object<-load_fatemap_into_seurat(data_dir, doublets_pct=act_dbl)
    prefix1 <- basename(data_dir)
    prefix2 <- tail(strsplit(data_dir, "/")[[1]], 3) |>
      head(1)
    doublet_object_file <- glue("{save_dir}/{prefix2}___{prefix1}___exp_{exp_dbl}__act_{act_dbl}.rds")
    cur_res<-run_doublet_finder(seurat_object, exp_dbl, doublet_object_file)
    df_ls<-rbindlist(list(df_ls,cur_res))
  }
  all_stats<-data.frame(df_ls)
  all_stats<-cbind(ID=rep(data_sample_ID,nrow(all_stats)), all_stats)
  all_stats
}

run_doublet_finder_on_dataset<-function(data_dir, METHOD_ID, exp_dbl, actual_dbl, num_to_test=5,
                                        save_root = "/projects/b1042/GoyalLab/zzhang/doublet_objects/doublet_finder"){
  print(glue("Working on dataset: {data_dir}"))
  dir_10X <- glue("{data_dir}/10X")
  dataset_id <- strsplit(data_dir,"/")[[1]][[length(strsplit(data_dir,"/")[[1]])]]
  sample_dirs<-list.dirs(path = dir_10X, full.names = TRUE, recursive = F)
  stats <- list()
  out_df<-NULL
  for (cur_sample_dir in sample_dirs){
    cur_sample_id<-strsplit(cur_sample_dir,"/")[[1]][[length(strsplit(cur_sample_dir,"/")[[1]])]]
    cur_data_sample_id<-glue("{dataset_id}_{cur_sample_id}")
    cur_save_dir <- glue("{save_root}/act__{actual_dbl}")
    if (!file.exists(cur_save_dir)) {
      dir.create(cur_save_dir, recursive = TRUE)
    }
    cur_df <- run_bchmk(cur_sample_dir, cur_data_sample_id, num_to_test, exp_dbl = exp_dbl,
                        act_dbl = actual_dbl, save_dir = cur_save_dir)
    out_df<-rbindlist(list(out_df, cur_df))|>
      data.frame()
  }
  out_df
}
