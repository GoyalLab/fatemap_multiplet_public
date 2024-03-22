
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
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final_sum"
save_root<-"/projects/b1042/GoyalLab/zzhang/doublet_objects_sum/doublet_finder"
num_test <- 1
doublet_rates_to_test<-c(0.05, 0.08, 0.1, 0.15, 0.2, 0.25)
ACT_DBL<-0.08
################################
#                              #
#             Main             #
#                              #
################################

all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
out_df<-NULL

for(cur_dataset_dir in all_datasets_dirs){
  print(cur_dataset_dir)
  # print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{ACT_DBL}"))
  out1 <- run_doublet_finder_on_dataset(cur_dataset_dir, "doublet_finder", actual_dbl = ACT_DBL, exp_dbl = 0.05, num_to_test = 1, save_root = save_root)
  out1$dbl_exp <- 0.05
  out1$dbl_act <- ACT_DBL
  out2 <- run_doublet_finder_on_dataset(cur_dataset_dir, "doublet_finder", actual_dbl = ACT_DBL, exp_dbl = 0.08, num_to_test = 1, save_root = save_root)
  out2$dbl_exp <- 0.08
  out2$dbl_act <- ACT_DBL
  out3 <- run_doublet_finder_on_dataset(cur_dataset_dir, "doublet_finder", actual_dbl = ACT_DBL, exp_dbl = 0.1, num_to_test = 1, save_root = save_root)
  out3$dbl_exp <- 0.1
  out3$dbl_act <- ACT_DBL
  out4 <- run_doublet_finder_on_dataset(cur_dataset_dir, "doublet_finder", actual_dbl = ACT_DBL, exp_dbl = 0.15, num_to_test = 1, save_root = save_root)
  out4$dbl_exp <- 0.15
  out4$dbl_act <- ACT_DBL
  out5 <- run_doublet_finder_on_dataset(cur_dataset_dir, "doublet_finder", actual_dbl = ACT_DBL, exp_dbl = 0.2, num_to_test = 1, save_root = save_root)
  out5$dbl_exp <- 0.2
  out5$dbl_act <- ACT_DBL
  out6 <- run_doublet_finder_on_dataset(cur_dataset_dir, "doublet_finder", actual_dbl = ACT_DBL, exp_dbl = 0.25, num_to_test = 1, save_root = save_root)
  out6$dbl_exp <- 0.25
  out6$dbl_act <- ACT_DBL

  cur_out_df <- rbindlist(list(out1, out2, out3, out4, out5, out6))|>
    data.frame()
  out_df<-rbindlist(list(out_df, cur_out_df))|>
    data.frame()
  print(out_df)
}

all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_{ACT_DBL}.tsv")
dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)