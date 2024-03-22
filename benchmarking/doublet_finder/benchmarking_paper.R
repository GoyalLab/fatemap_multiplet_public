
require(DoubletFinder)
require(PRROC)
require(Matrix)
require(Seurat)
require(glue)
require(purrr)
require(data.table)
require(dplyr)
require(stringr)
set.seed(2022)
source("doublet_finder_utils.R")

################################
#                              #
#       Input Parameters       #
#                              #
################################


root_dir <-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
root_bchmk_dataset<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
METHOD_ID<-"doublet_finder"
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final"
num_test <- 1
################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir, n_cores, num_test){
  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL

  print("INFO: Running on benchmark datasets!")

  all_bchmk_dataset_mtx<-list.files(path=root_bchmk_dataset, pattern="*.0.08.mtx", recursive = T, full.names = T)
  for(cur_bchmk_dataset_mtx in all_bchmk_dataset_mtx){
    df_ls<-NULL
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
      cur_dbl_rate <- sum(doublet.seurat[["label"]] == "doublet")/ncol(doublet.seurat)
      print(cur_dbl_rate)
      cur_file <- glue("/projects/b1042/GoyalLab/zzhang/doublet_objects_benchmarking_data/doublet_finder/{cur_data_id}.rds")
      cur_res<-run_doublet_finder(doublet.seurat, exp_dbl = cur_dbl_rate, doublet_obj_file = cur_file)
      df_ls<-rbindlist(list(df_ls,cur_res))
    }
    all_stats<-data.frame(df_ls)
    cur_out_df<-cbind(ID=rep(cur_data_id,nrow(all_stats)), all_stats)
    cur_out_df$dbl_exp <- "NA"
    cur_out_df$dbl_act <- "NA"
    out_df<-rbindlist(list(out_df, cur_out_df))|>
      data.frame()
    print("INFO: Done with original paper datasets")

  }


  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_original_benchmarking.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir, n_cores, num_test)
}