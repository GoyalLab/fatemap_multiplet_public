set.seed(2022)
library(Matrix)
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

source("scDblFinder_utils.R")
options(error=traceback)
################################
#                              #
#       Input Parameters       #
#                              #
################################


root_dir <-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
root_bchmk_dataset<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
METHOD_ID<-"scDblFinder"
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final"
n_cores<-8
num_test <- 1

################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir, n_cores){

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
  main(root_dir, METHOD_ID, out_dir, n_cores)
}
