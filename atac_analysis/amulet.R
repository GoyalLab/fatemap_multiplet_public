library(scDblFinder)

fragments <- "/projects/b1042/GoyalLab/zzhang/atac_data/cellTag/atac_fragments.tsv.gz"
res <- amulet(fragments)
write.table(res, file = "/projects/b1042/GoyalLab/zzhang/atac_data/new.tsv", quote = FALSE, sep="\t",row.names = TRUE)
