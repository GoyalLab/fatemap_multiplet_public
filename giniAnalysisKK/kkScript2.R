source(util.R)
root <- "/Users/karunkiani/Desktop/doublet/"
setwd(root)

loadLibs()


theme_set(theme_classic())
set.seed(123)

dataDir <- here('data')
plotDir <- here('plots')

Gini_tib <- tibble(
  Gini_Coeff = numeric(), 
  sample = character(),
  resolution = character()
)  

## Load in from 10X output files and create Seurat objects
####FM01####

tmpDir <- here(dataDir, 'FM01', '10X')

##FM01 - sample 3
tmp <- Read10X(data.dir = here(tmpDir, 'sample3'))
FM01_sample3 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM01_sample3"

barcodes <- read_tsv(here(dataDir, 'FM01', 'barcode', 'FM01_singlets.txt'),
                     col_names = FALSE)

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
  
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:50)
  seurat.obj <- FindClusters(seurat.obj, resolution = c(0.4,0.8,1.2),
                             save.SNN =TRUE)
  
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:50)
  
  
  Gini_tib <- makePlots(seurat.obj, sample_name, Gini_tib)
  
  return(Gini_tib) 
  
}

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

## Singlets function matches the seurat cellIDs with a list of barcodes we read in that we know to be from singlets
singlets <- function(seurat.obj, barcodes){
  
  seurat.obj@meta.data <- seurat.obj@meta.data %>% 
    dplyr::mutate(cellID = rownames(seurat.obj@meta.data)) %>% 
    tidyr::separate(cellID, into = c("cellID", NA)) %>% 
    dplyr::mutate(singlet = case_when(cellID %in% barcodes$X1 ~ "singlet",
                                      TRUE ~ "not singlet"))
  return(seurat.obj)
}  

Gini_tib <- runSample(FM01_sample3, barcodes, sample_name, Gini_tib)

##FM01 - sample 4
tmp <- Read10X(data.dir = here(tmpDir, 'sample4'))
FM01_sample4 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM01_sample4"

Gini_tib <- runSample(FM01_sample4, barcodes, sample_name, Gini_tib)

####FM02####

tmpDir <- here(dataDir, 'FM02', '10X')

barcodes <- read_tsv(here(dataDir, 'FM02', 'barcode', 'FM02_singlets.txt'),
                     col_names = FALSE)

##1uMPLX
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix'))
FM02_1uMPLX <- CreateSeuratObject(counts = tmp)
sample_name <- "FM02_1uMPLX"

Gini_tib <- runSample(FM02_1uMPLX, barcodes, sample_name, Gini_tib)

##100nMPLX
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix2'))
FM02_100nMPLX <- CreateSeuratObject(counts = tmp)
sample_name <- "FM02_100nMPLX"

Gini_tib <- runSample(FM02_100nMPLX, barcodes, sample_name, Gini_tib)

##5nMT
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix3'))
FM02_5nMT <- CreateSeuratObject(counts = tmp)
sample_name <- "FM02_5nMT"

Gini_tib <- runSample(FM02_5nMT, barcodes, sample_name, Gini_tib)

##5nMT100nMP
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix4'))
FM02_5nMT100nMP <- CreateSeuratObject(counts = tmp)
sample_name <- "FM02_5nMT100nMP"

Gini_tib <- runSample(FM02_5nMT100nMP, barcodes, sample_name, Gini_tib)


####FM03####

tmpDir <- here(dataDir, 'FM03', '10X')
barcodes <- read_tsv(here(dataDir, 'FM03', 'barcode', 'FM03_singlets.txt'),
                     col_names = FALSE)

##DMSO_1uMPLX
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix'))
FM03_DMSO_1uMPLX <- CreateSeuratObject(counts = tmp)
sample_name <- "FM03_DMSO_1uMPLX"

Gini_tib <- runSample(FM03_DMSO_1uMPLX, barcodes, sample_name, Gini_tib)

##DOT1Li_1uMPLX
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix2'))
FM03_DOT1Li_1uMPLX <- CreateSeuratObject(counts = tmp)
sample_name <- "FM03_DOT1Li_1uMPLX"

Gini_tib <- runSample(FM03_DOT1Li_1uMPLX, barcodes, sample_name, Gini_tib)


####FM04####

tmpDir <- here(dataDir, 'FM04', '10X')
barcodes <- read_tsv(here(dataDir, 'FM04', 'barcode', 'FM04_singlets.txt'),
                     col_names = FALSE)

##BC18_B1
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix'))
FM04_BC18_B1 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM04_BC18_B1"

Gini_tib <- runSample(FM04_BC18_B1, barcodes, sample_name, Gini_tib)

##BC18_B2
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix2'))
FM04_BC18_B2 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM04_BC18_B2"

Gini_tib <- runSample(FM04_BC18_B2, barcodes, sample_name, Gini_tib)

####FM05####

tmpDir <- here(dataDir, 'FM05', '10X')
barcodes <- read_tsv(here(dataDir, 'FM05', 'barcode', 'FM05_singlets.txt'),
                     col_names = FALSE)

##250nMPLX-A1
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix'))
FM05_250nMPLX_A1 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM05_250nMPLX_A1"

Gini_tib <- runSample(FM05_250nMPLX_A1, barcodes, sample_name, Gini_tib)

##250nMPLX-A2
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix2'))
FM05_250nMPLX_A2 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM05_250nMPLX_A2"

Gini_tib <- runSample(FM05_250nMPLX_A2, barcodes, sample_name, Gini_tib)

####FM06####

tmpDir <- here(dataDir, 'FM06', '10X')
barcodes <- read_tsv(here(dataDir, 'FM06', 'barcode', 'FM06_singlets.txt'),
                     col_names = FALSE)

##WM989Naive_1
tmp <- Read10X(data.dir = here(tmpDir,  'FM06-WM989Naive-1', 'filtered_feature_bc_matrix'))
FM06_WM989Naive_1 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM06_WM989Naive_1"

Gini_tib <- runSample(FM06_WM989Naive_1, barcodes, sample_name, Gini_tib)

##WM989Naive_2
tmp <- Read10X(data.dir = here(tmpDir, 'FM06-WM989Naive-2', 'filtered_feature_bc_matrix'))
FM06_WM989Naive_2 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM06_WM989Naive_2"

Gini_tib <- runSample(FM06_WM989Naive_2, barcodes, sample_name, Gini_tib)

####FMXX####

tmpDir <- here(dataDir, 'FMXX', '10X')
barcodes <- read_tsv(here(dataDir, 'FMXX', 'barcode', 'FMXX_singlets.txt'),
                     col_names = FALSE)

##sample3
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix'))
FMXX_sample3 <- CreateSeuratObject(counts = tmp)
sample_name <- "FMXX_sample3"

Gini_tib <- runSample(FMXX_sample3, barcodes, sample_name, Gini_tib)
##sample4
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix2'))
FMXX_sample4 <- CreateSeuratObject(counts = tmp)
sample_name <- "FMXX_sample4"

Gini_tib <- runSample(FMXX_sample4, barcodes, sample_name, Gini_tib)

##sample5
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix3'))
FMXX_sample5 <- CreateSeuratObject(counts = tmp)
sample_name <- "FMXX_sample5"

Gini_tib <- runSample(FMXX_sample5, barcodes, sample_name, Gini_tib)

####NJ01####

tmpDir <- here(dataDir, 'NJ01', '10X')
barcodes <- read_tsv(here(dataDir, 'NJ01', 'barcode', 'NJ01_singlets.txt'),
                     col_names = FALSE)

##DMSO_A
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix'))
NJ01_DMSO_A <- CreateSeuratObject(counts = tmp)
sample_name <- "NJ01_DMSO_A"

Gini_tib <- runSample(NJ01_DMSO_A, barcodes, sample_name, Gini_tib)

##DMSO_B
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix2'))
NJ01_DMSO_B <- CreateSeuratObject(counts = tmp)
sample_name <- "NJ01_DMSO_B"

Gini_tib <- runSample(NJ01_DMSO_B, barcodes, sample_name, Gini_tib)

##LSD1i_A
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix3'))
NJ01_LSD1i_A <- CreateSeuratObject(counts = tmp)
sample_name <- "NJ01_LSD1i_A"

Gini_tib <- runSample(NJ01_LSD1i_A, barcodes, sample_name, Gini_tib)

##LSDi_B
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix4'))
NJ01_LSD1i_B <- CreateSeuratObject(counts = tmp)
sample_name <- "NJ01_LSD1i_B"

Gini_tib <- runSample(NJ01_LSD1i_B, barcodes, sample_name, Gini_tib)

##DOT1Li_A
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix5'))
NJ01_DOT1Li_A <- CreateSeuratObject(counts = tmp)
sample_name <- "NJ01_DOT1Li_A"

Gini_tib <- runSample(NJ01_DOT1Li_A, barcodes, sample_name, Gini_tib)

##DOT1Li_B
tmp <- Read10X(data.dir = here(tmpDir, 'filtered_feature_bc_matrix6'))
NJ01_DOT1Li_B <- CreateSeuratObject(counts = tmp)
sample_name <- "NJ01_DOT1Li_B"

Gini_tib <- runSample(NJ01_DOT1Li_B, barcodes, sample_name, Gini_tib)

####Make Distribution Plot for Gini tibble####

pal <- park_palette("Saguaro", 4)

plot <- ggboxplot(Gini_tib, y = "Gini_Coeff", x = "resolution",
                  color = "resolution", palette = c(pal[2:4]),
                  add = "jitter", shape = "resolution") +
  rotate()

p <- gghistogram(Gini_tib,  x = "Gini_Coeff",
                 add='mean', rug=TRUE,
                 color="resolution", fill ="resolution",
                 palette = pal[1:3], facet.by ="resolution",alpha = .9)