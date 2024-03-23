library(dplyr)
library(MASS)
library(ggrepel)
library(matrixStats)

library(ggplot2)

find_peak <- function(data) {
  # Compute the kernel density estimate on a grid
  kde <- MASS::kde2d(data$umap_1, data$umap_2, n = 1000)
  
  # Find the grid point with the highest density
  max_density_index <- which(kde$z == max(kde$z), arr.ind = TRUE)
  
  # Return the coordinates of the highest density
  peak_coords <- cbind(kde$x[max_density_index[,1]], kde$y[max_density_index[,2]])
  
  # Return a data frame with peak coordinates
  return(data.frame(umap_1 = peak_coords[1, 1], umap_2 = peak_coords[1, 2]))
}

# Function to read a CSV file and plot its data, saving the plot to a specified directory
plot_csv <- function(file_path, plot_dir, plot_png) {
  # Ensure the plot directory exists
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Read the CSV file
  df <- read.csv(file_path)
  peaks <- do.call(rbind, by(df, df$cluster, find_peak))
  peaks$cluster = rownames(peaks)
  # Assuming the CSV has columns named 'x' and 'y'
  umap = ggplot(df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, shape = 16) +
    geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
                     box.padding = unit(1, "lines"), point.padding = unit(1, "lines"),
                      force = 5, family="Arial",max.overlaps = Inf, 
                     color = "black", bg.color = "white", bg.r = 0.08, size = 6,min.segment.length	= 0.05,nudge_y = -1, nudge_x = 1)  +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umap)
  
  #Constructing path for the plot
    plot_file_name <- paste0(tools::file_path_sans_ext(basename(file_path)), ".svg")
    plot_file_path <- file.path(plot_dir, plot_file_name)
  # Save the plot
    svglite(filename = plot_file_path, width = 6, height = 4)
    plot(umap)
    dev.off()
    
    png_file_name <- paste0(tools::file_path_sans_ext(basename(file_path)), ".png")
    png_file_path <- file.path(plot_png, png_file_name)
    ggsave(umap, file = png_file_path, width = 6, height = 4)
}

# Function to read a CSV file and plot its data, saving the plot to a specified directory
replot_csv <- function(file_path, plot_dir, plot_png) {
  # Ensure the plot directory exists
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Read the CSV file
  df <- read.csv(file_path)
  peaks <- do.call(rbind, by(df, df$cluster, find_peak))
  peaks$cluster = rownames(peaks)
  # Assuming the CSV has columns named 'x' and 'y'
  umap = ggplot(df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, shape = 16) +
    geom_text_repel(data = peaks, aes(x = umap_1, y = umap_2, label = cluster), 
                    box.padding = unit(1, "lines"), point.padding = unit(1, "lines"),
                    force = 100, family="Arial",max.overlaps = Inf, max.time = 20,
                    color = "black", bg.color = "white", bg.r = 0.08, size = 6,min.segment.length	= 0.05)  +
    theme_void() +
    theme(legend.position = "none", 
          plot.background = element_blank(), # This sets the plot background to be blank
          panel.background = element_blank(), # Ensure the panel background is blank
          panel.border = element_blank(), # Removes panel border if present
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank()) # Removes minor grid lines
  
  # To display the plot (if running interactively)
  print(umap)
  
  #Constructing path for the plot
  plot_file_name <- paste0(tools::file_path_sans_ext(basename(file_path)), ".svg")
  plot_file_path <- file.path(plot_dir, plot_file_name)
  # Save the plot
  svglite(filename = plot_file_path, width = 6, height = 4)
  plot(umap)
  dev.off()
  
  png_file_name <- paste0(tools::file_path_sans_ext(basename(file_path)), ".png")
  png_file_path <- file.path(plot_png, png_file_name)
  ggsave(umap, file = png_file_path, width = 6, height = 4)
}

# Recursive function to iterate through all files in a directory and its subdirectories
plot_all_csv <- function(dir_path, plot_dir, plot_png) {
  # List all files and directories at the current path
  contents <- list.files(dir_path, full.names = TRUE)
  
  # Iterate over each item
  for (item in contents) {
    dataset = tools::file_path_sans_ext(basename(dir_path))
    if (dir.exists(item)) {
      # If the item is a directory, recurse into it
      plot_all_csv(item, plot_dir, plot_png)
    } else if (grepl("\\.csv$", item)) {
      # If the item is a CSV file, plot it and save the plot to the specified directory
      print(paste("Plotting and saving:", item))
      plot_csv(item, plot_dir, plot_png)
    }
  }
}

# Base directory to start searching for CSV files and target directory for plots
base_dir <- "~/Keerthana/CellTypeCountData/plotData/UmapData" # Update this path
plot_dir <- "~/Keerthana/CellTypeCountData/plots/PlotsFormatted" # Update this path
plot_png <- "~/Keerthana/CellTypeCountData/plots/pngPlots"

plot_all_csv(base_dir, plot_dir, plot_png)
rePlotList <- c("/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D3-RNA-r1-1.csv", 
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D3-RNA-r1-2.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D3-RNA-r2-4.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D3-RNA-r2-5.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D21-RNA-r1-1.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D21-RNA-r1-2.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D21-RNA-r1-4.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/CellTag/CellTag_B4D21-RNA-r2-3.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/ClonMapper/ClonMapper_FM1.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/ClonMapper/ClonMapper_FM7.csv",
                "/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/LARRY/LARRY_LK1_d4_R_4.csv"
                )
for (csvFile in rePlotList){
  print(paste("Plotting and saving:", csvFile))
  replot_csv(csvFile,plot_dir,plot_png)
}
plot_csv("/home/mzo5929/Keerthana/CellTypeCountData/plotData/UmapData/FateMap-NonCancer",plot_dir,plot_png)
