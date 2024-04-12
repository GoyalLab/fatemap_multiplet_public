require(Seurat)
require(glue)
require(purrr)
require(data.table)
require(dplyr)
require(stringr)



seed <- 2022
set.seed(seed)

s_root <- "/projects/p31666/zzhang/doublet-bchmk/data/functional_analysis_dataset"
# output_root <- "/projects/p31666/zzhang/doublet-bchmk/data/DE_output_supp_run"
output_root <- "/projects/p31666/zzhang/doublet-bchmk/data/DE_output"
# TODO: Change if rerunning all
# all_dataset_dir <- list.dirs(s_root, recursive = F)
all_dataset_dir <- c(
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/Biorxiv",
  # "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/LARRY",
  "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/ClonMapper"
  # "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/SPLINTR",
  # "/projects/b1042/GoyalLab/zzhang/functional_dataset_overflow/TREX"
)
for(cur_dataset_dir in all_dataset_dir){
  cur_dataset_id <- basename(cur_dataset_dir)
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
  Idents(forty_seu) <- "seurat_clusters"
  cur_clusters <- Idents(forty_seu)|>levels()
  combinations <- combn(cur_clusters, 2)
  num_combinations_to_sample <- min(3, ncol(combinations))
  selected_combinations <- combinations[, sample(ncol(combinations), num_combinations_to_sample)]

  cur_dataset_out_dir <- glue("{output_root}/{seed}/{cur_dataset_id}")
  dir.create(cur_dataset_out_dir, recursive = TRUE, showWarnings = FALSE)
  for(cur_dataset_dbl_pct in names(cur_dataset_all_s)){
    cur_s <- cur_dataset_all_s[[cur_dataset_dbl_pct]]
    cur_s[["RNA"]] <- JoinLayers(cur_s[["RNA"]])
    Idents(cur_s) <- "seurat_clusters"
    for(cur_col in ncol(selected_combinations)){
      cur_clus_1 <- selected_combinations[1, cur_col]
      cur_clus_2 <- selected_combinations[2, cur_col]

      cur_out_file <- glue("{cur_dataset_out_dir}/{cur_dataset_dbl_pct}_{cur_clus_1}_vs._{cur_clus_2}__wilcox.csv")
      if(file.exists(cur_out_file) == FALSE){
        cur_DE <- FindMarkers(cur_s, ident.1 = cur_clus_1, ident.2 = cur_clus_2)
        cur_DE[["gene"]] <- rownames(cur_DE)
        write.csv(cur_DE, file = cur_out_file, quote = F, row.names = F)
      }else{
        print(glue("Skip: {cur_out_file}"))
      }

      cur_out_file <- glue("{cur_dataset_out_dir}/{cur_dataset_dbl_pct}_{cur_clus_1}_vs._{cur_clus_2}__MAST.csv")
      if(file.exists(cur_out_file) == FALSE){
        cur_DE <- FindMarkers(cur_s, ident.1 = cur_clus_1, ident.2 = cur_clus_2, test.use = "MAST")
        cur_DE[["gene"]] <- rownames(cur_DE)
        write.csv(cur_DE, file = cur_out_file, quote = F, row.names = F)
      }else{
        print(glue("Skip: {cur_out_file}"))
      }
    }
  }
}
