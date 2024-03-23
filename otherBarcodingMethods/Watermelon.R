#Importing the required libraries
library(Seurat)
#BiocManager::install("ShortRead")
library("ShortRead")
library(DropletUtils)
library(data.table)
library(stringr)


#Reading in the Seurat file for all samples
T47D <- readRDS(file = "~/Keerthana/Watermelon_Atac-seq/T47D_seurat.rds")
genes <- Features(T47D)

valid.10Xcell.barcodes.df=read.table("/home/mzo5929/Keerthana/Watermelon_Atac-seq/Watermelon_otherData/10X/barcodes.tsv.gz",sep = "-")
validBarcodes10X=valid.10Xcell.barcodes.df$V1

#create cell names as metadata colum
T47D[["CellName"]] <- Cells(T47D)
#output.path = "~/Keerthana/Watermelon_Atac-seq/Watermelon_otherData/10X/"
#write10xCounts(x = T47D@assays$RNA@layers$counts, path = output.path, gene.id = genes, barcodes = cells)

###############################################################################
#Helper function to split the seurat object based on the sample (using Cell IDs) and save it in a separate folder in 10X format (matrix.mtx, genes.tsv.gz, barcodes.tsv.gz)
splitSaveSeuratObject <- function(sampleName, CellName){
  cat("Subsetting the Seurat object\n")
  CellName <- paste0(CellName,"-1")
  cat("Number of final cells", length( CellName), "\n")
  seurat_subset <- subset(T47D, subset = Cells(T47D) %in% CellName)
  

  cat("Creating a separate directory for the sample if it already not present\n")
  sample_folder <- file.path(output_path, sampleName)
  dir.create(sample_folder, showWarnings = TRUE)
  
  cat("Save Seurat object in the sample folder\n")
  write10xCounts(x = seurat_subset@assays$RNA@layers$counts, path = sample_folder, gene.id = genes, barcodes = Cells(seurat_subset))
  
}
###############################################################################
CollapseLineageBarcodes <- function(single_10x_bc.df){
  if (nrow(single_10x_bc.df) == 0) {
    cat("Input data frame is empty.\n")
    return(single_10x_bc.df)
  }
  #how many umi support this 10x barcode
  original_umi_count=sum(single_10x_bc.df$count)
  #order such that the most supported pair would be the first
  single_10x_bc.df=single_10x_bc.df[order(-single_10x_bc.df$count),]
  #create a vector of all lineage barcodes
  str<- single_10x_bc.df$barcode
  #Levenshtein Distance of all lineage barcodes
  str<- single_10x_bc.df$barcode
  d  <- adist(str)
  #add a col with edit distance to most common lineage barcode
  single_10x_bc.df$edist_to_most_common_LB= d[1, ]
  #go over data frame and remove any row that has an edit distance of 3 or less
  single_10x_bc.df=single_10x_bc.df[!(single_10x_bc.df$edist_to_most_common_LB<4 & single_10x_bc.df$edist_to_most_common_LB>0),]
  #calculate how many umi were subtracted and add these umis to the most common barcode
  single_10x_bc.df$count[1]=single_10x_bc.df$count[1]+original_umi_count-sum(single_10x_bc.df$count)
  
  #if after collapsing the less supported lineages are supported only by one barcode and the leading barcode is larger than 5-remove these rows
  #mark top count
  top_count=single_10x_bc.df$count[1]
  
  #remove all pairs where the support in 3 times smaller compared to the top pair
  single_10x_bc.df=single_10x_bc.df[!(single_10x_bc.df$count*3<top_count),]  
  return(single_10x_bc.df)
}

###############################################################################
#Using the fastq files to recreate Barcode mapping similar to lineage_db.csv and saving the seurat object for each sample
Barcode_10X_mapping<-function(sampleName, read1, read2, validBarcodes10X){
  #Open the fastq file
  R1_fq_da <-readFastq(dir_path, pattern=read1) #dial out R1 read
  R2_fq_da <-readFastq(dir_path, pattern=read2) #dial out R2 read
  
  #Extract the sequences - R1, R2, cell Id
  #extracting sequences
    #Extracting read 1 and read 2 from respective fastq files
    R1_seq_da=as.character(R1_fq_da@sread) #for dial-out
    R2_seq_da=as.character(R2_fq_da@sread) #for dial-out
    
    #creating a dataframe to store the information from miseq data
    dialOut.df=as.data.frame(cbind(R1_seq_da,R2_seq_da)) 
    
    #extracting cell ID from read 1
    dialOut.df$cellID =gsub("^([ATCG]{16}).*","\\1",dialOut.df$R1,perl=TRUE )

    #Making the cell ID format similar to the one in a Seurat object by adding sample ID(typically 1)
    #dialOut.df$cellID=paste0(dialOut.df$cellID,"-1")
  cat("Just from reads", nrow(dialOut.df),"\n")  
  
  #Match the cell barcodes found in this to the ones in 10X data
  #matching cell ID in this data to 10X
  cat(length(unique(dialOut.df$cellID)), "\n")
  dialOut.df=subset(dialOut.df, dialOut.df$cellID %in% valid.10Xcell.barcodes.df$V1)
  
  
  #Merge multiple reads of same (Cell ID, Lineage Barcode, UMI) - make each row unique essentially
  #Merging multiple read
  dialOut.df=unique(dialOut.df)
  cat("num unique", nrow(dialOut.df), "\n")
  
  #Extract UMI, Lineage Barcode  from the read sequences
  #extracting UMI and lineage barcode
    #extracting umi of each read (last 10bp of read R1)
  dialOut.df$umi=gsub("^[ATCG]{16}(.*)","\\1",dialOut.df$R1,perl=TRUE )
    
    #Extracting the lineage barcode based on the pattern it is supposed to have
    pattern <- "GGGCTG(([AT][CG]|[CG][AT]){15})GACGCT"
    sequences <- as.character(dialOut.df$R2_seq_da)
    matches <- str_match(dialOut.df$R2_seq_da, pattern)
    dialOut.df$barcode<- matches[,2]
  
  #Filter based on valid lineage barcodes
    dialOut.df=subset(dialOut.df,nchar(dialOut.df$barcode) == 30)
  cat("Valid Barcodes", nrow(dialOut.df),"\n")  

  
  #Calculate the number of UMIs for every pair of Cell Barcode - Lineage Barcode
  #Calculating UMIs
  dialOut_2=setDT(dialOut.df)[, .(count = uniqueN(umi)), by =list(cellID,barcode) ]
  cat("calculating UMIs ", nrow(dialOut_2), "\n")
  #Collapse the barcodes using Levenshtein distance and other parameters(which barcode to merge to and what should be the maximum edit distance to merge, etc) using the pre-defined helper function
  #Collapsing based on LV distance
  dialOut_2<-setDT(dialOut_2)[, CollapseLineageBarcodes(.SD), by=cellID, .SDcols=c("barcode", "count")]
  cat("Final collapse", nrow(dialOut_2),"\n") 
  
  
  #Subset the Seurat object using the cell IDs and save it in the necessary 10X format in a separate folder using a predefined function
  #Saving it in 10X format
  splitSaveSeuratObject(sampleName, dialOut_2$cellID)
  
  #Return the list of Cell ID, Lineage Barcode, UMI counts and Sample so that it can be merged in the full list
  dialOut_2$sample <- sampleName
  finalData <- dialOut_2[, c("cellID", "barcode", "count", "sample"), drop = FALSE]
  cat("final data:", sampleName,dim(dialOut_2), "\n")
  return(finalData) 
}
#####################
#The dataframe to store cellID, Lineage Barcode, UMI counts and samples for all the samples
cellID_LB_UMI_Sample <- data.frame(
  cellID = character(),
  barcode = character(),
  UMICount = numeric(),
  sample = character(),
  stringsAsFactors = FALSE
)

# FASTQ files to be input
filenames_fastq <- c(
  "T47D-lag-2_S4_L001_R1_001.fastq.gz", "T47D-late-2_S6_L001_R1_001.fastq.gz", "T47D-naive-2_S2_L001_R1_001.fastq.gz",
  "T47D-lag-2_S4_L001_R2_001.fastq.gz", "T47D-late-2_S6_L001_R2_001.fastq.gz", "T47D-naive-2_S2_L001_R2_001.fastq.gz",
  "T47D-lag-1_S3_L001_R1_001.fastq.gz", "T47D-late-1_S5_L001_R1_001.fastq.gz", "T47D-naive-1_S1_L001_R1_001.fastq.gz",
  "T47D-lag-1_S3_L001_R2_001.fastq.gz", "T47D-late-1_S5_L001_R2_001.fastq.gz", "T47D-naive-1_S1_L001_R2_001.fastq.gz"
)

dir_path="~/Keerthana/Watermelon_Atac-seq/Watermelon_otherData/Fastq/" #fastq dir path for dial out 
output_path="~/Keerthana/Watermelon_Atac-seq/Watermelon_otherData/Processed_Fastq/"     #path for final output file

#filenames_fastq <- c("T47D-lag-1_S3_L001_R1_001.fastq.gz", "T47D-lag-1_S3_L001_R2_001.fastq.gz")

# Sort filenames to ensure pairs are adjacent
sorted_filenames <- sort(filenames_fastq)

# Process two filenames at a time
for (i in seq(1, length(sorted_filenames), by = 2)) {
  # Extract the two filenames
  filename_1 <- sorted_filenames[i]
  filename_2 <- sorted_filenames[i + 1]

  sample_name <- sapply(strsplit(filename_1, "_"), function(x) x[1])
  
  # Print or process the pair as needed
  cat("Processing sample:", sample_name, "\n")
  
  # Add your processing logic here
  temp_data <- Barcode_10X_mapping(sample_name, filename_1, filename_2, validBarcodes10X)
  cellID_LB_UMI_Sample <- rbind(cellID_LB_UMI_Sample, temp_data)
}

df_expanded <- cellID_LB_UMI_Sample[rep(seq_len(nrow(cellID_LB_UMI_Sample)), cellID_LB_UMI_Sample$count), ]
sum(cellID_LB_UMI_Sample$count)
# Drop the 'count' column
df_expanded <- df_expanded[, -3, drop = FALSE]

# Assuming 'df' is your dataframe
write.csv(df_expanded, file = paste0(output_path,"Watermelon_barcode_cellID_UMIcounts.csv"), row.names = FALSE)









