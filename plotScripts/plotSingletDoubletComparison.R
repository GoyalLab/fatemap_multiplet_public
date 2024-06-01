library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

newCompiledNumbers <- read_csv("Keerthana/CellCounts/newCompiledNumbers.csv")
newCompiledNumbers$ExperimentName <- gsub("FateMap_FM(\\d+)", "Goyal et al.\\1", newCompiledNumbers$ExperimentName)

#Plotting number of barcodes for samples after QC vs singletCode recovery after QC
barcodedCell <- newCompiledNumbers$numBarcodedCellsQC
singletCodeCell <- newCompiledNumbers$numTrueSingletsQC

singletCodevsBarcodedLegend <- ggplot(newCompiledNumbers, aes(x = numBarcodedCellsQC, y = numTrueSingletsQC, color = ExperimentName)) +
  geom_point(size = 2, shape = 16) +
  scale_color_viridis_d(option = "viridis") + # Discrete viridis color scale
  theme_minimal() +
  theme(legend.position = "right", 
        plot.background = element_blank(), # This sets the plot background to be blank
        panel.background = element_blank(),
        panel.grid.minor = element_blank()) # Ensure the panel background is blank


ggsave("~/Keerthana/CellCounts/singletCodevsBarcoded.svg", plot = singletCodevsBarcoded, height = 4, width = 6)        
ggsave("~/Keerthana/CellCounts/singletCodevsBarcodedwithLegend.png", plot = singletCodevsBarcodedLegend, height = 4, width = 6) 

singletCodevsBarcoded <- ggplot(newCompiledNumbers, aes(x = numBarcodedCellsQC, y = numTrueSingletsQC, color = ExperimentName)) +
  geom_point(size = 2, shape = 16) +
  # scale_color_viridis_d(option = "viridis") + # Discrete viridis color scale
  theme_minimal() +
  theme(legend.position = "none", 
        plot.background = element_blank(), # This sets the plot background to be blank
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_blank() ) # Ensure the panel background is blank

singletCodevsBarcoded
ggsave("~/Keerthana/CellCounts/singletCodevsBarcoded.svg", plot = singletCodevsBarcoded, height = 4, width = 4)        
ggsave("~/Keerthana/CellCounts/singletCodevsBarcoded.png", plot = singletCodevsBarcoded, height = 4, width = 4) 


#Plotting comparable singlet recovery rate 
experiments <- c("Goyal et al.01", "Goyal et al.02", "Goyal et al.03","Goyal et al.04",
                 "Goyal et al.05", "Goyal et al.06", "Goyal et al.08", "JainEtAl", "JiangEtAl", "SPLINTR")
comparableExperiments <- newCompiledNumbers %>% filter(ExperimentName %in% experiments)
data_long <- pivot_longer(comparableExperiments, cols = c("numPaperSingletsQC", "numTrueSingletsQC"), names_to = "MeasurementType", values_to = "Value")
average_data <- data_long %>%
  group_by(MeasurementType) %>%
  summarise(AverageValue = mean(Value))%>% ungroup()

singletCodevsPaper <- ggplot() +
  # Plot points for each sample with color by experiment, in the background
  geom_point(data = data_long, aes(x = MeasurementType, y = Value, color = ExperimentName), 
             size = 1, shape = 16, alpha = 0.5) +
  # Connect points within each sample with gray lines, also in the background
  geom_line(data = data_long, aes(x = MeasurementType, y = Value, group = sampleNameList), 
            alpha = 0.1, color = "gray") +
  # Add average line, bolder and on top
  geom_point(data = average_data, aes(x = MeasurementType, y = AverageValue, color = "black"), 
             size = 2, shape = 16) +
  geom_line(data = average_data, aes(x = MeasurementType, y = AverageValue, group = 1), 
            linewidth = 1, color = "black") + # Adjust color and size as desired for the average line
  theme_minimal() +
  labs() +
  theme(legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

singletCodevsPaper
ggsave("~/Keerthana/CellCounts/singletCodevsPaperAveraged.svg", plot = singletCodevsPaper, height = 4, width = 4)        
ggsave("~/Keerthana/CellCounts/singletCodevsPaperAveraged.png", plot = singletCodevsPaper, height = 4, width = 4) 

singletCodevsPaperLegend <- ggplot() +
  # Plot points for each sample with color by experiment, in the background
  geom_point(data = data_long, aes(x = MeasurementType, y = Value, color = ExperimentName), 
             size = 1, shape = 16, alpha = 0.5) +
  # Connect points within each sample with gray lines, also in the background
  geom_line(data = data_long, aes(x = MeasurementType, y = Value, group = sampleNameList), 
            alpha = 0.1, color = "gray") +
  # Add average line, bolder and on top
  geom_point(data = average_data, aes(x = MeasurementType, y = AverageValue, color = "black"), 
             size = 2, shape = 16) +
  geom_line(data = average_data, aes(x = MeasurementType, y = AverageValue, group = 1), 
            linewidth = 1, color = "black") + # Adjust color and size as desired for the average line
  theme_minimal() +
  labs() 
singletCodevsPaperLegend
ggsave("~/Keerthana/CellCounts/singletCodevsPaperLegendAveraged.svg", plot = singletCodevsPaperLegend, height = 4, width = 6)        
ggsave("~/Keerthana/CellCounts/singletCodevsPaperLegendAveraged.png", plot = singletCodevsPaperLegend, height = 4, width = 6) 

paperSingletAvg <- mean(comparableExperiments$numPaperSingletsQC)
singletCodeSingletAvg <- mean(comparableExperiments$numTrueSingletsQC)

singletCodevsPaperLegend <- ggplot(data_long, aes(x = MeasurementType, y = Value, group = sampleNameList, color = ExperimentName)) +
  geom_point(size = 2, shape = 16) +
  #scale_color_manual(values = color) +
  geom_line() +
  labs() +
  theme_minimal() 
singletCodevsPaperLegend
ggsave("~/Keerthana/CellCounts/singletCodevsPaperLegend.svg", plot = singletCodevsPaperLegend, height = 4, width = 6)        
ggsave("~/Keerthana/CellCounts/singletCodevsPaperLegend.png", plot = singletCodevsPaperLegend, height = 4, width = 6) 
