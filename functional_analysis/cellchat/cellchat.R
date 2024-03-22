library(CellChat)
require(Seurat)
require(glue)
require(purrr)
require(data.table)
require(dplyr)
require(stringr)


s_root <- "/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset"
s_overflow_root <- "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow"
output_root <- "/projects/p31666/zzhang/doublet-bchmk/data/cellchat_output"
# all_dataset_dir <- list.dirs(s_root, recursive = F)
# partial_rerun
all_dataset_dir <- c(
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/Biorxiv",
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/LARRY",
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/ClonMapper",
  # "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/SPLINTR",
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/TREX"
)
for(cur_dataset_dir in all_dataset_dir){
  print(cur_dataset_dir)
  cur_dataset_id <- basename(cur_dataset_dir)

  if(grepl("smart", cur_dataset_id)){
    print(glue("Manually skipped {cur_dataset_id}"))
    next
  }

  if(cur_dataset_id %in% c("cellTag", "TREX", "LARRY")){
    CellChatDB <- CellChatDB.mouse
  }else{
    CellChatDB <- CellChatDB.human
  }

  # overflow case
  if(cur_dataset_id %in% c("TREX", "SPLINTR", "Biorxiv", "LARRY", "ClonMapper")){
    cur_dataset_dir <- str_replace(cur_dataset_dir, s_root, s_overflow_root)
  }
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

  if(file.exists(glue("{cur_dataset_dir}/data/forty.rds")) == FALSE){
    next
  }
  forty_seu <- readRDS(glue("{cur_dataset_dir}/data/forty.rds"))
  twenty_seu <- readRDS(glue("{cur_dataset_dir}/data/twenty.rds"))
  ten_seu <- readRDS(glue("{cur_dataset_dir}/data/ten.rds"))
  control_seu <- readRDS(glue("{cur_dataset_dir}/data/control.rds"))
  cur_dataset_all_s <- list(
    "forty" = forty_seu,
    "twenty" = twenty_seu,
    "ten" = ten_seu,
    "control" = control_seu
  )

  cur_dataset_out_dir <- glue("{output_root}/{cur_dataset_id}")
  dir.create(cur_dataset_out_dir, recursive = TRUE, showWarnings = FALSE)

  for(cur_dataset_dbl_pct in names(cur_dataset_all_s)){
    cur_out_file <- glue("{cur_dataset_out_dir}/{cur_dataset_dbl_pct}.csv")
    if(file.exists(cur_out_file)){
      print(glue("{cur_out_file} exists. Skip to next!"))
      next
    }
    cur_s <- cur_dataset_all_s[[cur_dataset_dbl_pct]]
    cur_s[["RNA"]] <- JoinLayers(cur_s[["RNA"]])
    cur_s@meta.data[["clusters_c"]] <- paste0("C", cur_s$seurat_clusters)
    Idents(cur_s) <- "clusters_c"
    cellchat <- createCellChat(object = cur_s)
    if(grepl("SPLINTR", cur_dataset_id)){
      mtx <- cellchat@data
      rownames(mtx) <- toupper(rownames(mtx))
      cellchat@data <- mtx
    }
    cellchat@DB <- CellChatDB.use
    cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    df.net <- subsetCommunication(cellchat)

    write.csv(df.net, file = cur_out_file, quote = F, row.names = F)
  }
}