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
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final_sum"
save_root<-"/projects/b1042/GoyalLab/zzhang/doublet_objects_sum/scDblFinder"
n_cores<-8
doublet_rates_to_test<-c(0.05, 0.08, 0.1, 0.15, 0.2, 0.25)
num_test <- 1
ACT_DBL <- 0.05

################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir, n_cores){

  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL

  # loop through doublet rates to test for
    for(expected_doublet_rate in doublet_rates_to_test){
      for(cur_dataset_dir in all_datasets_dirs){
        print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{ACT_DBL}"))
        cur_out_df<-run_scDblFinder_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
                                                dbl_exp = expected_doublet_rate, dbl_act = ACT_DBL, save_root = save_root)

        cur_out_df$dbl_exp <- expected_doublet_rate
        cur_out_df$dbl_act <- ACT_DBL
        print(cur_out_df)
        out_df<-rbindlist(list(out_df, cur_out_df))|>
          data.frame()
      }
    }

  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_{ACT_DBL}.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir, n_cores)
}
