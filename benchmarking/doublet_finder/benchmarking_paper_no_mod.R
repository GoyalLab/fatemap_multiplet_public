
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
  out_df<-NULL

  print("INFO: Running on benchmark datasets!")

  all_bchmk_dataset<-list.files(path=root_bchmk_dataset, pattern="*.rds", recursive = T, full.names = T)
  for(cur_bchmk_dataset in all_bchmk_dataset){
    df_ls<-NULL
    for(idx in 1:num_test){
      cur_data<-readRDS(cur_bchmk_dataset)
      cur_data_id<-strsplit(cur_bchmk_dataset, "/")[[1]][[8]]
      count<-cur_data[[1]]
      label<-cur_data[[2]]
      doublet.seurat <- CreateSeuratObject(counts = count, min.cells = 1)
      doublet.seurat <- AddMetaData(doublet.seurat, col.name="label", label)
      cur_dbl_rate <- sum(doublet.seurat[["label"]] == "doublet")/ncol(doublet.seurat)
      cur_file <- glue("/projects/b1042/GoyalLab/zzhang/doublet_objects_benchmarking_data_no_mod/doublet_finder/{cur_data_id}.rds")
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


  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_original_benchmarking_no_mod.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir, n_cores, num_test)
}