data_dir<-"/projects/p31666/zzhang/doublet-bchmk/data/couting_threshold/scale_factor_0.00003"
out_dir <- "/projects/p31666/zzhang/doublet-bchmk/plots/dataset_vis"
all_out_files<-list.files(data_dir, pattern = "*singlets_stats.csv", full.names = T, recursive = T)
all_two_sample_barcode_singlets_files <- list.files(data_dir, pattern = "*all_sample_two_barcode_singlets.txt", full.names = T, recursive = T)
# all_cells <- sum(all_stats_df$total_cells)
# all_singlets <- sum(all_stats_df$total_singlets)
# print(all_singlets / all_cells)

all_stats_df <- NULL

for(cur_file in all_out_files){
  cur_df <- read.csv(cur_file)
  all_stats_df <- rbindlist(list(all_stats_df, cur_df))
}
all_stats_df <- all_stats_df|>as.data.frame()
all_stats_df_singlets <- dplyr::select(all_stats_df,
                                       dataset,
                                       single_sample_barcode_singlets,
                                       multi_barcode_singlets,
                                       dominant_umi_barcode_singlets)

all_stats_df_long <- gather(all_stats_df_singlets, category, count, -dataset)

# add two sample singlets
cur_stats_df <- list()
for(cur_file in all_two_sample_barcode_singlets_files){
  cur_dataset <- sub(".*/(.*?)_all.*", "\\1", cur_file)
  print(cur_dataset)
  if (file.info(cur_file)$size == 0) {
    cur_count <- 0
  }else{
    cur_df <- read.csv(cur_file, header = F)
    cur_count <- nrow(cur_df)
    print(length(cur_df$V1) == length(unique(cur_df$V1)))
  }
  cur_stats_df[[cur_file]] <- list(
    "dataset" = cur_dataset,
    "category" = "multi_sample_barcode_singlets",
    "count" = cur_count
  )
}
cur_stats_df <- rbindlist(cur_stats_df)|>
  as.data.frame()
all_stats_df_long_combined <- rbind(all_stats_df_long, cur_stats_df)
all_stats_df_long_combined <- dplyr::filter(all_stats_df_long_combined, dataset != "smartseq3")

svg(glue("{out_dir}/singlets_by_category.svg"), width = 15, height = 8)
p<-ggplot(all_stats_df_long_combined, aes(x = dataset, y = count, fill = category)) +
  geom_bar(stat = "identity")+
  theme_Publication_blank()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(p)
dev.off()

all_stats_df_long_pct <- all_stats_df_long_combined %>%
  dplyr::group_by(dataset, category) %>%
  dplyr::summarize(count = sum(count), .groups = 'drop')%>%
  group_by(dataset) %>%
  dplyr::mutate(percentage = count / sum(count) * 100)

p <- ggplot(all_stats_df_long_pct, aes(x = percentage, y = dataset, fill = category)) +
  geom_bar(stat = "identity") +
  theme_Publication_blank() +
  ylab("Percentage") # Label the y-axis as "Percentage"
svg(glue("{out_dir}/singlets_by_category_pct.svg"), width = 15, height = 8)
plot(p)
dev.off()
write.csv(all_stats_df_long_pct, glue("{out_dir}/singlets_by_category_pct.csv"), quote = F, row.names = F)