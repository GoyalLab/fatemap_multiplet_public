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
source("scDblFinder_utils.R")
################################
#                              #
#       Input Parameters       #
#                              #
################################


root_dir <-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
METHOD_ID<-"scDblFinders"
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final"
root_bchmk_dataset<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir){
  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-list()

  all_bchmk_dataset<-list.files(path=root_bchmk_dataset, pattern="*.rds", recursive = T, full.names = T)
  for(cur_bchmk_dataset in all_bchmk_dataset){
    cur_data<-readRDS(cur_bchmk_dataset)
    cur_data_id<-strsplit(cur_bchmk_dataset, "/")[[1]][[8]]
    print(cur_data_id)
    count<-cur_data[[1]]
    label<-cur_data[[2]]
    doublet.seurat <- CreateSeuratObject(counts = count, min.cells = 1)
    doublet.seurat <- AddMetaData(doublet.seurat, col.name="label", label)
    sce <- as.SingleCellExperiment(doublet.seurat)
    cur_exp_dbl<-sum(label=="doublet")/length(label)
    cur_file <- glue("/projects/b1042/GoyalLab/zzhang/doublet_objects_benchmarking_data_no_mod/scDblFinders/{cur_data_id}.rds")
    cur_out_df<-run_scDblFinder(sce, dbl_exp=cur_exp_dbl, doublet_obj_file=cur_file)
    print(glue("INFO: Current benchmark dataset doublet rate is: {cur_exp_dbl}"))
    cur_out_df$dbl_exp <- "NA"
    cur_out_df$dbl_act <- "NA"
    cur_out_df$ID <- cur_data_id
    out_df[[cur_data_id]] <- cur_out_df
  }
  out_df<-rbindlist(out_df)|>
    as.data.frame()

  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_original_benchmarking_no_mod.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir)
}
