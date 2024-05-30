loadLibs <- function(){
  library(tidyverse)
  library(Seurat)
  library(scCustomize)
  library(DescTools)
  library(nationalparkcolors)
  library(ggpubr)
  library(ggridges)
  library(here)
  
}


#Takes a seurat object and barcodes dataframe and outputs the seurat object
#with an extra metadata column that identifies which cellIDs are in the list
#of singlet barcodes
singlets <- function(seurat.obj, barcodes){
  
  seurat.obj@meta.data <- seurat.obj@meta.data %>% 
    dplyr::mutate(cellID = rownames(seurat.obj@meta.data)) %>% 
    tidyr::separate(cellID, into = c("cellID", NA)) %>% 
    dplyr::mutate(singlet = case_when(cellID %in% barcodes$X1 ~ "singlet",
                                      TRUE ~ "not singlet"))
  return(seurat.obj)
}  



#takes the seurat object, sample_name (character), data frame of Gini Coeffs
#plots QC metrics, seurat clusters, and plots of which cells are singlets for multiple resolutions
#also plots the proportion of singlets per seurat cluster and calculates Gini Coefficients of proportions
#and adds the relevant data to the Gini tibble
makePlots <- function(seurat.obj, sample_name, Gini_tib){
  
  resolutions <- c("RNA_snn_res.0.4",
                   "RNA_snn_res.0.8",
                   "RNA_snn_res.1.2")
  
  pal <- park_palette("Arches", 4)
  
  for (i in 1:3){
    
    res <- resolutions[i]
    color <- pal[i+1]
    
    Idents(seurat.obj) <- res
    
    #make violin plot of QC data
    
    VlnPlot_scCustom(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                     group.by = res)
    
    ggsave(file= here(plotDir,paste0(sample_name,"_",res,"_violinQC.svg")))
    
    ##make umap of singlets
    
    DimPlot_scCustom(seurat.obj, group.by = "singlet",
                     colors_use = c(pal[1], color),
                     pt.size = 0.75) +
      ggtitle(paste0(res, " ", "singlets"))
    
    ggsave(file= here(plotDir,paste0(sample_name,"_",res,"_singlets.svg")))
    
    ##make umap of clusters
    
    DimPlot_scCustom(seurat.obj, reduction = 'umap', 
                     label = TRUE,
                     pt.size = 0.75) +
      ggtitle(paste0(res, " ", "clusters"))
    
    ggsave(file = here(plotDir,paste0(sample_name,"_",res,"_seuratClust.svg")))
    
    ##Make proportion plot
    plotTib <- as.data.frame(table(Idents(seurat.obj), seurat.obj@meta.data$singlet))
    plotTib <- as_tibble(plotTib) 
    
    plotTib <- plotTib %>% 
      dplyr::rename(cluster = Var1) %>% 
      dplyr::rename(singlet = Var2)
    
    plot <- ggplot(plotTib, aes(fill=singlet, y=Freq, x=cluster)) +
      geom_bar(position = "fill", stat = "identity") + 
      scale_fill_manual(values = c(pal[1], color)) +
      xlab("Seurat cluster") +
      ylab("Proportion of cells") + 
      theme(legend.position = "top") +
      theme(legend.title=element_blank()) +
      ggtitle(paste0(sample_name, " ", res, " ", "singlet proportions by cluster"))
    
    ggsave(file = here(plotDir,paste0(sample_name,"_",res,"_proportionBarPlot.svg")))
    
    ##calculate Gini Coeff and add result to table
    
    tmp_tib <- seurat.obj@meta.data %>% 
      group_by(!!as.symbol(res), singlet) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% 
      filter(singlet == "singlet")
    
    freq <- tmp_tib$freq
    
    Gini_tib <- Gini_tib %>%
      add_row(Gini_Coeff = Gini(freq), 
              sample = sample_name, 
              resolution = res) 
  }
  
  return(Gini_tib)
  
}

#Wrapper Script that takes a seurat object, barcodes data frame, sample name (char), and Empty 
#Gini tibble which will be used to graph results and performs QC and filtering of seurat object
#as well as running PCA and UMAP and uses singlets and makePlots helper functions
runSample <- function(seurat.obj, barcodes, sample_name, Gini_tib){
  
  #Run singlets function
  seurat.obj <- singlets(seurat.obj, barcodes)
  
  ## Preprocessing steps
  
  #mito genes
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  
  feature_quartile <- quantile(seurat.obj@meta.data$nFeature_RNA,
                               probs = seq(0.25, 0.75, by = 0.25))
  
  mt_decile <- quantile(seurat.obj@meta.data$percent.mt,
                        probs = seq(0.1, 0.9, by = 0.1))
  
  #filter
  seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > feature_quartile[1]
                       & percent.mt < mt_decile[6])
  
  #Normalize
  seurat.obj <- NormalizeData(seurat.obj)
  
  #Highly variable features
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  
  #scale data on variable features
  all.genes <- rownames(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj, features = all.genes)
  
  #Dimension reduction
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  
  #Write PC data to dataframe
  out <- Embeddings(seurat.obj, reduction = "pca")
  
  #add cellID to data frame
  out <- cbind(rownames(out), as.data.frame(out, row.names = NULL))
  colnames(out)[1] <- "cellID"
  
  #write to TSV
  write_tsv(x = out, file = paste0(here("extractedData", "PC_tables", ""), sample_name, "_PCtable.tsv"))
  
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:50)
  seurat.obj <- FindClusters(seurat.obj, resolution = c(0.4,0.8,1.2),
                             save.SNN =TRUE)
  
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:50)
  
  
  Gini_tib <- makePlots(seurat.obj, sample_name, Gini_tib)
  
  return(Gini_tib) 
  
}
