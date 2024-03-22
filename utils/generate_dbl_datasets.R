require(Seurat)
require(glue)
require(purrr)
require(data.table)
require(dplyr)
require(stringr)

root_path<-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
# out_root <- "/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset"
out_root <- "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow"

# all_directories <- list.dirs(root_path, full.names = TRUE, recursive = TRUE)
# all_directories <-  grep("10X$", all_directories, value = TRUE)
all_directories <- c(
  "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/SPLINTR/10X",
  "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/TREX/10X"
)

load_fatemap_into_seurat <- function(data_dir, doublets_pct = 0.08) {
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
  singlet_count <- sum(doublet.seurat@meta.data$label == "singlet")
  # Seurat V5 modifications
  t <- JoinLayers(doublet.seurat)
  doublet.seurat[["RNA"]] <- t[["RNA"]]
  doublet.seurat
}



for(cur_10X_dir in all_directories){
  print(cur_10X_dir)
  all_10X_dir <- list.dirs(cur_10X_dir, full.names = T, recursive = F)
  forth_seu_ls <- list()
  parts <- unlist(strsplit(cur_10X_dir, "/"))
  dataset_id <- parts[length(parts) - 1]
  cur_dataset_plot_dir <- glue("{out_root}/{dataset_id}/plots")
  cur_dataset_data_dir <- glue("{out_root}/{dataset_id}/data")
  dir.create(cur_dataset_plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cur_dataset_data_dir, recursive = TRUE, showWarnings = FALSE)

  for(cur_10X_sample_dir in all_10X_dir){
    add_id <- basename(cur_10X_sample_dir)

    cur_data_40 <- load_fatemap_into_seurat(cur_10X_sample_dir, doublets_pct = 0.4)
    cur_data_40 <- RenameCells(cur_data_40, add.cell.id = add_id)
    cur_data_40[["sample"]] <- add_id
    forth_seu_ls[[cur_10X_sample_dir]] <- cur_data_40
  }

  forty_seu <- merge(forth_seu_ls[[1]], y = forth_seu_ls[-1])
  forty_seu[["RNA"]] <- JoinLayers(forty_seu[["RNA"]])
  forty_seu <- NormalizeData(forty_seu)
  forty_seu <- ScaleData(forty_seu)
  forty_seu <- FindVariableFeatures(forty_seu, selection.method = "vst", nfeatures = 2000)
  forty_seu <- RunPCA(forty_seu)
  forty_seu[["RNA"]] <- split(forty_seu[["RNA"]], f = forty_seu$sample)
  forty_seu <- IntegrateLayers(object = forty_seu, method = CCAIntegration,
                               orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
  forty_seu <- FindNeighbors(forty_seu, reduction = "integrated.cca", dims = 1:30)
  forty_seu <- FindClusters(forty_seu, resolution = 0.5)
  forty_seu <- RunUMAP(forty_seu, dims = 1:30, reduction = "integrated.cca")


  p1 <- DimPlot(forty_seu, raster = F)
  forty_cluster_pdf <- glue("{cur_dataset_plot_dir}/forty_cluster.pdf")
  pdf(forty_cluster_pdf)
  plot(p1)
  dev.off()


  Idents(forty_seu) <- "label"
  p2 <- DimPlot(forty_seu, raster = F)
  forty_label_pdf <- glue("{cur_dataset_plot_dir}/forty_label.pdf")
  pdf(forty_label_pdf)
  plot(p2)
  dev.off()

  Idents(forty_seu) <- "sample"
  p3 <- DimPlot(forty_seu, raster = F)
  forty_sample_pdf <- glue("{cur_dataset_plot_dir}/forty_sample.pdf")
  pdf(forty_sample_pdf)
  plot(p3)
  dev.off()
  forty_data_file <- glue("{cur_dataset_data_dir}/forty.rds")
  saveRDS(forty_seu, file = forty_data_file)

  doublet_cells <- forty_seu@meta.data|>
    dplyr::filter(label == "doublet")|>
    rownames()
  singlet_cells <- forty_seu@meta.data|>
    dplyr::filter(label == "singlet")|>
    rownames()
  num_doublets <- round(length(singlet_cells) * 0.2 / (1 - 0.2))
  #randomly sample without replacement
  doublet_to_keep <- sample(doublet_cells, num_doublets)
  all_cells_twenty <- c(singlet_cells, doublet_to_keep)
  twenty_seu <- forty_seu[,all_cells_twenty]

  Idents(twenty_seu) <- "seurat_clusters"
  p1 <- DimPlot(twenty_seu, raster = F)
  twenty_cluster_pdf <- glue("{cur_dataset_plot_dir}/twenty_cluster.pdf")
  pdf(twenty_cluster_pdf)
  plot(p1)
  dev.off()


  Idents(twenty_seu) <- "label"
  p2 <- DimPlot(twenty_seu, raster = F)
  twenty_label_pdf <- glue("{cur_dataset_plot_dir}/twenty_label.pdf")
  pdf(twenty_label_pdf)
  plot(p2)
  dev.off()

  num_doublets <- round(length(singlet_cells) * 0.1 / (1 - 0.1))
  #randomly sample without replacement
  doublet_to_keep <- sample(doublet_cells, num_doublets)
  all_cells_ten <- c(singlet_cells, doublet_to_keep)
  ten_seu <- forty_seu[,all_cells_ten]

  Idents(ten_seu) <- "seurat_clusters"
  p1 <- DimPlot(ten_seu, raster = F)
  ten_cluster_pdf <- glue("{cur_dataset_plot_dir}/ten_cluster.pdf")
  pdf(ten_cluster_pdf)
  plot(p1)
  dev.off()


  Idents(ten_seu) <- "label"
  p2 <- DimPlot(ten_seu, raster = F)
  ten_label_pdf <- glue("{cur_dataset_plot_dir}/ten_label.pdf")
  pdf(ten_label_pdf)
  plot(p2)
  dev.off()



  control_seu <- subset(forty_seu, label == "singlet")
  Idents(control_seu) <- "seurat_clusters"
  p1 <- DimPlot(control_seu, raster = F)
  control_cluster_pdf <- glue("{cur_dataset_plot_dir}/control_cluster.pdf")
  pdf(control_cluster_pdf)
  plot(p1)
  dev.off()

  Idents(control_seu) <- "label"
  p2 <- DimPlot(control_seu, raster = F)
  control_label_pdf <- glue("{cur_dataset_plot_dir}/control_label.pdf")
  pdf(control_label_pdf)
  plot(p2)
  dev.off()

  twenty_data_file <- glue("{cur_dataset_data_dir}/twenty.rds")
  saveRDS(twenty_seu, file = twenty_data_file)

  ten_data_file <- glue("{cur_dataset_data_dir}/ten.rds")
  saveRDS(ten_seu, file = ten_data_file)

  control_data_file <- glue("{cur_dataset_data_dir}/control.rds")
  saveRDS(control_seu, file = control_data_file)
}
