#set root directory where you will access raw files and output tables and plots
root <- "/Users/karunkiani/Desktop/doublet/"
setwd(root)

#point to the other script where the functions are kept
source(paste0(root, "extractionScripts/util.R"))

#Load necessary librarires, set themes, set seed
loadLibs()
theme_set(theme_classic())
set.seed(123)

#set data and plot directories
dataDir <- here('data')
plotDir <- here('plots')

#Create blank Gini Tibble to start off with
Gini_tib <- tibble(
  Gini_Coeff = numeric(), 
  sample = character(),
  resolution = character()
)  

### Work through each sample to run the wrapper function ####
  
## Load in from 10X output files and create Seurat objects
####FM01####

tmpDir <- here(dataDir, 'FM01', '10X')

barcodes <- read_tsv(here(dataDir, 'FM01', 'barcode', 'FM01_singlets.txt'),
                     col_names = FALSE)


##FM01 - sample 3
tmp <- Read10X(data.dir = here(tmpDir, 'sample3'))
FM01_sample3 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM01_sample3"

Gini_tib <- runSample(FM01_sample3, barcodes, sample_name, Gini_tib)

##FM01 - sample 4
tmp <- Read10X(data.dir = here(tmpDir, 'sample4'))
FM01_sample4 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM01_sample4"

Gini_tib <- runSample(FM01_sample4, barcodes, sample_name, Gini_tib)

##FM01 - sample 1
tmp <- Read10X(data.dir = here(tmpDir, 'sample1'))
FM01_sample1 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM01_sample1"

Gini_tib <- runSample(FM01_sample1, barcodes, sample_name, Gini_tib)

##FM01 - sample 2
tmp <- Read10X(data.dir = here(tmpDir, 'sample2'))
FM01_sample2 <- CreateSeuratObject(counts = tmp)
sample_name <- "FM01_sample2"

Gini_tib <- runSample(FM01_sample2, barcodes, sample_name, Gini_tib)



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

####CJ01####

tmpDir <- here(dataDir, 'CJ01', '10X')
barcodes <- read_tsv(here(dataDir, 'CJ01', 'barcode', 'CJ01_singlets.txt'),
                     col_names = FALSE)

##Barsibling_A
tmp <- Read10X(data.dir = here(tmpDir, '10kbarsiblingA'))
CJ01_Barsibling_A <- CreateSeuratObject(counts = tmp)
sample_name <- "CJ01_Barsibling_A"

Gini_tib <- runSample(CJ01_Barsibling_A, barcodes, sample_name, Gini_tib)

##Barsibling_B
tmp <- Read10X(data.dir = here(tmpDir, '10kbarsiblingB'))
CJ01_Barsibling_B <- CreateSeuratObject(counts = tmp)
sample_name <- "CJ01_Barsibling_B"

Gini_tib <- runSample(CJ01_Barsibling_B, barcodes, sample_name, Gini_tib)



### Write out Gini-tib ####
write_tsv(x = Gini_tib, file = paste0(here("extractedData",""), "Gini_Plot_data.tsv"))

#########################################


####Make Distribution Plot for Gini tibble####

p <- gghistogram(Gini_tib,  x = "Gini_Coeff",
            add='mean', rug=TRUE,
            color="resolution", fill ="resolution",
            palette = pal[1:3], facet.by ="resolution",alpha = .9) +
            xlim(c(0,1))

ggsave(p, file = "/Users/karunkiani/Desktop/doublet/plots/giniDistributions.svg",
       dpi = 300, height = 6, width = 12, units = "in")



### Create histograms of example distributions and their respective ####
## Gini distributiosn


#create empty distribution tibble
Dist_tib <- tibble(
  Gini_Coeff = numeric(), 
  sample = character(),
  distribution = character()
)  

#Function to find Gini distribution for multiple iterations of different dists
Gini_dist <- function(Dist_tib){
  for (i in 1:23){
    distribution <- "uniform"
    nums <- runif(n = 20, min = 0, max = 1)
  
    sample = paste0(distribution, "_", i)
    
    Dist_tib <- Dist_tib %>%
      add_row(Gini_Coeff = Gini(nums), 
              sample = paste0(distribution, "_", i), 
              distribution = distribution) 
  }
  
  for (i in 1:23){
    distribution <- "exponential"
    nums <- rexp(n = 20)
    
    sample = paste0(distribution, "_", i)
    
    Dist_tib <- Dist_tib %>%
      add_row(Gini_Coeff = Gini(nums), 
              sample = paste0(distribution, "_", i), 
              distribution = distribution) 
  }
  
  for (i in 1:23){
    distribution <- "power"
    nums <- poweRlaw::rplcon(20,xmin = 1, 1.5)
    
    sample = paste0(distribution, "_", i)
    
    Dist_tib <- Dist_tib %>%
      add_row(Gini_Coeff = Gini(nums), 
              sample = paste0(distribution, "_", i), 
              distribution = distribution) 
  }

return(Dist_tib)

}

Dist_tib <- Gini_dist(Dist_tib)

Dist_tib <- Dist_tib %>% 
  mutate(distribution = fct_relevel(distribution,
                                    "uniform",
                                    "exponential",
                                    "power"))

p2 <- gghistogram(Dist_tib,  x = "Gini_Coeff",
                 add='mean', rug=TRUE,
                 color="distribution", fill ="distribution",
                 facet.by ="distribution",alpha = .9) + xlim(c(0,1))

ggsave(p2, file = "/Users/karunkiani/Desktop/doublet/plots/exampleDistributions.svg",
       dpi = 300, height = 6, width = 12, units = "in")

dat <- as.data.frame(runif(n = 10000, min = 0, max = 1))
colnames(dat) <- "val"

p3 <- gghistogram(dat, x = "val",
                  add = "mean", rug = TRUE, 
                  color = "dodgerblue", fill = "dodgerblue",
                  alpha=0.9)

ggsave(p3, file = "/Users/karunkiani/Desktop/doublet/plots/uniform.svg",
       dpi = 300)

dat <- as.data.frame(rexp(n =100000))
colnames(dat) <- "val"

p4 <- gghistogram(dat, x = "val",
                  add = "mean", rug = TRUE, 
                  color = "dodgerblue", fill = "dodgerblue",
                  alpha=0.9)

ggsave(p4, file = "/Users/karunkiani/Desktop/doublet/plots/exponential.svg",
       dpi = 300)

dat <- as.data.frame(poweRlaw::rplcon(10000,xmin = 1, 1.5))
colnames(dat) <- "val"

p5 <- gghistogram(dat, x = "val",
                  add = "mean", rug = TRUE, 
                  color = "dodgerblue", fill = "dodgerblue",
                  alpha=0.9) + xlim(c(0, 100))
ggsave(p5, file = "/Users/karunkiani/Desktop/doublet/plots/power.svg",
       dpi = 300)
