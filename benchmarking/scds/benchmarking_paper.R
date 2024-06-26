library(Matrix)
library(Seurat)
library(scds)
library(SingleCellExperiment)
library(PRROC)
require(glue)
require(data.table)
library(stringr)
require(parallel)

set.seed(2022)
source("scds_utils.R")
################################
#                              #
#       Input Parameters       #
#                              #
################################


root_dir <-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
METHOD_ID<-"hybrid"
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final"
root_bchmk_dataset<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
num_test <- 1
ACT_DBL <- 0.05
################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir){
  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL

  all_bchmk_dataset_mtx<-list.files(path=root_bchmk_dataset, pattern="*.0.08.mtx", recursive = T, full.names = T)
  for(cur_bchmk_dataset_mtx in all_bchmk_dataset_mtx){
    cur_out_df<-run_bchmk_original_paper(cur_bchmk_dataset_mtx, num_to_test = num_test)
    cur_out_df$dbl_exp <- "NA"
    cur_out_df$dbl_act <- "NA"
    out_df<-rbindlist(list(out_df, cur_out_df))|>
      data.frame()
  }

  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_original_benchmarking.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir)
}
