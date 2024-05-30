#Steps: 
#read in the curve points and separate it by lineages
#read in the linear trajectory
#do similarity matching and calculate penalty score
#after matching, find lock-step euclidean distance and multiply by penalty
#normalise for total number of lineages
install.packages("ggpubr")
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(clue)
library(ggpubr)
# Function Definitions
calculateEditDistance <- function(path1, path2) {
  # Ensure paths are numeric
  path1 <- as.numeric(path1)
  path2 <- as.numeric(path2)
  
  # Initialize matrix
  len1 <- length(path1)
  len2 <- length(path2)
  matrix <- matrix(0, nrow = len1 + 1, ncol = len2 + 1)
  
  # Set up base conditions
  for (i in 1:(len1 + 1)) {
    matrix[i, 1] <- i - 1
  }
  for (j in 1:(len2 + 1)) {
    matrix[1, j] <- j - 1
  }
  
  # Compute distances
  for (i in 2:(len1 + 1)) {
    for (j in 2:(len2 + 1)) {
      cost <- ifelse(path1[i - 1] == path2[j - 1], 0, 1)
      matrix[i, j] <- min(matrix[i - 1, j] + 1, matrix[i, j - 1] + 1, matrix[i - 1, j - 1] + cost)
    }
  }
  
  # Return the Edit Distance
  return(matrix[len1 + 1, len2 + 1])
}

calculateCostMatrix <- function(controlLineages, sampleLineages) {
  # Assuming controlLineages and sampleLineages are lists of tibbles/data.frames as shown
  
  # Extracting numeric vectors from the tibble's Lineage colum
  
  # Initialize the cost matrix with dimensions [n x m]
  n <- length(controlLineages[[1]])
  m <- length(sampleLineages[[1]])
  costMatrix <- matrix(0, nrow = n, ncol = m)
  
  # Calculate pairwise edit distances to fill the cost matrix
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      # Using numeric vectors for edit distance calculation
      costMatrix[i, j] <- calculateEditDistance(controlLineages[[1]][[i]], sampleLineages[[1]][[j]])
    }
  }
  
  return(costMatrix)
}

matchLineages <- function(lineages, sampleMatch){
  control <- lineages[["control"]]
  match <- lineages[[sampleMatch]]
  costMatrix <- calculateCostMatrix(control, match)
  controlIsRow <- TRUE
  if (nrow(costMatrix) > ncol(costMatrix)) {
    costMatrix <- t(costMatrix)
    controlIsRow <- FALSE  # Indicates that control lineages are now columns, not rows
  }
  assignment <- solve_LSAP(costMatrix)
  len <- length(assignment)
  parsedMat <- cbind(
    1:len,
    as.integer(assignment)
  )
  if(controlIsRow) {
    colnames(parsedMat) <- c("Control", "Sample")
  } else {
    colnames(parsedMat) <- c("Sample", "Control")
  }
  parsedMat <- as.data.frame(parsedMat)
  # Correctly fetching edit distances based on whether control lineages were rows or columns
  if (controlIsRow) {
    parsedMat$editDistance <- mapply(function(i, j) costMatrix[i, j], parsedMat$Control, parsedMat$Sample)
  } else {
    parsedMat$editDistance <- mapply(function(i, j) costMatrix[j, i], parsedMat$Control, parsedMat$Sample)
  }
  costMatrix <- as.data.frame(costMatrix)
  distanceAndMatchedLineage <- list(costMatrix, parsedMat)
  return (distanceAndMatchedLineage)
}

# Function to process a single lineage row into a numeric vector
processLineageRow <- function(lineageRow) {
  as.numeric(unlist(strsplit(lineageRow, ",\\s*")))
}

#Not used
calculatePenalty <- function(diff) {
  # Simple absolute difference for penalty
  penalty <- abs(diff)
  return(penalty)
}

calculateLockstepEuclideanDistance <- function(trajectory1, trajectory2) {
  if (nrow(trajectory1) != nrow(trajectory2)) {
    stop("Trajectories are of different lengths!")
  }
  
  squaredDifferences <- rowSums((trajectory1 - trajectory2)^2)
  distance <- sqrt(sum(squaredDifferences))
  return(distance)
}

findDistanceAll <- function(experimentID){
  samples <- c("control", "ten", "twenty", "forty")
  inputTrajectoryFolder <- "~/Keerthana/CellTrajectory/data/curvesPlot/"
  inputLineageGraph <-  "~/Keerthana/CellTrajectory/lineageTrajectory/"
  
  # Initialize empty lists to store trajectories and lineages data
  trajectories <- list()
  lineages <- list()
  
  # Read in trajectory and lineage data
  for (sample in samples) {
    trajectories[[sample]] <- read_csv(paste0(inputTrajectoryFolder, experimentID, "_", sample, ".csv"))
    lineages[[sample]] <- read_csv(paste0(inputLineageGraph, experimentID, "_", sample, ".csv")) %>%
      mutate(Lineage = lapply(Lineage, function(x) as.numeric(unlist(strsplit(x, ",\\s*")))))
  }
  
  # Similarity Matching and Penalty Calculation
  matchedLineages <- list()
  editDistance <- list()
  for (sampleMatch in samples) { # Exclude 'control' from comparison
    distanceAndMatchedLineage <- matchLineages(lineages, sampleMatch)
    editDistance[[sampleMatch]] <- distanceAndMatchedLineage[[1]]
    matchedLineages[[sampleMatch]] <- distanceAndMatchedLineage[[2]]
  }
  for (sample in samples[-1]){
    write.csv(editDistance[[sample]], paste0("~/Keerthana/CellTrajectory/finalDistances/", experimentID,"_",sample, "_editDistance.csv"))
    
  }
  # Calculate Lockstep Euclidean Distance and Compile Results
  finalDistanceDf <- bind_rows(lapply(names(matchedLineages), function(sampleName) {
    mappingDf <- matchedLineages[[sampleName]]
    samplePairs <- list()
    
    for (i in 1:nrow(mappingDf)) {
      controlID <- mappingDf$Control[i]
      sampleID <- mappingDf$Sample[i]
      editDistance <- mappingDf$editDistance[i]
      controlTrajectory <- trajectories$control %>% filter(Lineage == controlID) %>% select(1:2,)
      numLineageControl <- length(unique(trajectories$control$Lineage))
      sampleTrajectory <- trajectories[[sampleName]] %>% filter(Lineage == sampleID) %>% select(1:2,)
      numLineageSample <- length(unique(trajectories[[sampleName]]$Lineage))
      distance <- calculateLockstepEuclideanDistance(controlTrajectory, sampleTrajectory)
      penalty <- calculatePenalty(length(unique(trajectories$control$Lineage)) - length(unique(trajectories[[sampleName]]$Lineage)))
      distanceWithPenalty <- distance*penalty
      samplePairs <- rbind(samplePairs, data.frame(SampleName = sampleName, numLineageControl = numLineageControl, numLineageSample = numLineageSample, ControlID = controlID, SampleID = sampleID, editDistance = editDistance, Distance = distance, Difference = penalty))
    }
    
    return(samplePairs)
  }), .id = "Sample")

  finalDistanceDf$ExperimentID <- experimentID
  return(finalDistanceDf)
}

################
#main
experimentList <- c("FM01", "FM02", "FM03", "FM04", "FM05", "FM06", "FM08", "Biorxiv", 
                    "cellTag", "ClonMapper","LARRY","non_cancer","smartseq3_umis","SPLINTR", "watermelon")
experimentResults <- lapply(experimentList, findDistanceAll)
masterTable <- bind_rows(experimentResults)
write.csv(masterTable, paste0("~/Keerthana/CellTrajectory/finalDistances/finalTable.csv"))
calculate_penalty <- function(difference) {
  # Simple example: a fixed penalty for any difference > 0
  if(difference != 0){
    penalty <-   exp(as.integer(difference))
    
  }else{
    penalty = 0
  }
  return(penalty)
}

comparisonTable <- masterTable %>%
  group_by(SampleName, ExperimentID) %>%
  summarise(
    AverageDistance = mean(Distance),  # Calculate average editDistance
    MaxDifference = max(Difference),  # Assuming you want the maximum difference in this group for penalty calculation
    Penalty = calculate_penalty(Difference[[1]])  # Apply penalty based on the max difference
  ) %>%
  ungroup() %>%
  mutate(
    AdjustedDistance = AverageDistance + Penalty  # Apply the penalty to the average distance
  )
comparisonTable <- comparisonTable %>%
  mutate(SampleName = factor(SampleName, levels = c("ten", "twenty", "forty")),  # Adjust if you have more or different samples
         ExperimentID = factor(ExperimentID, levels = unique(ExperimentID)))  # Preserves unique order or specify manually

library(ggplot2)
library(reshape2)


#Rough heatmap to try out
ggplot(comparisonTable, aes(x = ExperimentID, y = SampleName, fill = factor(AdjustedDistance))) +
  geom_tile(color = "white") +  # Use color to distinguish tiles
  scale_fill_viridis_d(begin = 0.3, end = 0.7, direction = -1, option = "A") +  # Discrete color scale
  labs(title = "Heatmap of Adjusted Distances",
       x = "Experiment ID",
       y = "Sample Name",
       fill = "Adjusted\nDistance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
#Function to account for penalties for differing number of lineages
calculate_penalty_other <- function(difference) {
  # Assuming you want to apply a penalty for any non-zero difference
  penalty <- ifelse(difference != 0, exp(difference), 0)  # Adjust the penalty calculation as needed
  return(penalty)
}

otherComparisonTable <- masterTable %>%
  mutate(
    Penalty = calculate_penalty_other(Difference),
    AdjustedDistance = Distance + 0.1*Distance*Penalty
  )
#Normalising the total distance to compare easily - max normalisation
otherComparisonTableNormalised <- otherComparisonTable %>%
  group_by(ExperimentID) %>%
  mutate(
    MaxAdjustedDistance = max(AdjustedDistance, na.rm = TRUE),  # Find the max AdjustedDistance within each group
    NormalizedAdjustedDistance = AdjustedDistance / MaxAdjustedDistance  # Normalize AdjustedDistance
  ) %>%
  ungroup()  %>%
  mutate(SampleName = factor(SampleName, levels = c("control","ten", "twenty", "forty")),  # Adjust if you have more or different samples
         ExperimentID = factor(ExperimentID, levels = unique(ExperimentID)))  # Preserves unique order or specify manually
#Normalisation usig max-min normalisation
otherComparisonTableNormalised <- otherComparisonTable %>%
  group_by(ExperimentID) %>%
  mutate(
    MinAdjustedDistance = min(AdjustedDistance, na.rm = TRUE),  # Find the min AdjustedDistance within each group
    MaxAdjustedDistance = max(AdjustedDistance, na.rm = TRUE),  # Find the max AdjustedDistance within each group
    NormalizedAdjustedDistance = (AdjustedDistance - MinAdjustedDistance) / (MaxAdjustedDistance - MinAdjustedDistance)  # Min-Max normalization
  ) %>%
  ungroup() %>%
  mutate(
    SampleName = factor(SampleName, levels = c("control","ten", "twenty", "forty")),  # Adjust if you have more or different samples
    ExperimentID = factor(ExperimentID, levels = unique(ExperimentID))  # Preserves unique order or specify manually
  )


write.csv(otherComparisonTableNormalised, paste0("~/Keerthana/CellTrajectory/finalDistances/finalTableWithPenalty.csv"))

otherComparisonTableNormalised = read_csv("~/doubletProject/tempTesting/finalTableWithPenalty.csv")

#Adding control rows to the plot (since it will match exactly with itself)
addControlRows <- function(df) {
  # Extract the numLineageControl value from the first row
  numLineageControl <- df$numLineageControl[1]
  
  # Create control rows
  controlRows <- df[rep(1, numLineageControl), ]
  controlRows$editDistance <- 0  # Set all columns after the Sample, SampleName, numLineageControl, and numLineageSample to zero
  controlRows$Distance  <- 0 
  controlRows$Difference <- 0 
  controlRows$Penalty <- 0 
  controlRows$AdjustedDistance <- 0 
  controlRows$NormalizedAdjustedDistance <- 0 
  controlRows$SampleID <- 1:numLineageControl  # Set sampleID from 1 to numLineageControl
  controlRows$ControlID <- 1:numLineageControl
  controlRows$numLineageSample <- controlRows$numLineageControl  # Set numLineageSample to numLineageControl
  controlRows$MaxAdjustedDistance <- df$MaxAdjustedDistance[1]
  controlRows$ExperimentID <- df$ExperimentID[1]
  controlRows$SampleName <- 'control'  # Set sampleName to 'control'
  controlRows$Sample <- 4
  controlRows$...1 <- df$...1[1]
  # Combine the original data with the control rows
  combinedData <- rbind(df, controlRows)
  
  return(combinedData)
}

# Apply the function to each sampleID and experimentID combination and combine
otherComparisonTableNormalised <- otherComparisonTableNormalised %>%
  group_by(ExperimentID) %>%
  do(addControlRows(.)) %>%
  ungroup()
otherComparisonTableNormalised$...1 <- rownames(otherComparisonTableNormalised)

median_data <- otherComparisonTableNormalised %>%
  group_by(ExperimentID, SampleName) %>%
  summarize(MedianNormalizedAdjustedDistance = median(NormalizedAdjustedDistance, na.rm = TRUE)) %>%
  ungroup()
# Plotting
color_blind_friendly_colors <- c("red", "lightblue", "steelblue", "darkblue")

allExperimentsBoxPlot <- ggplot(otherComparisonTableNormalised, aes(x = SampleName, y = NormalizedAdjustedDistance, fill = SampleName)) +
    geom_boxplot() +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
    facet_wrap(~ExperimentID, scales = "free_x") +
    geom_line(data = median_data, aes(x = SampleName, y = MedianNormalizedAdjustedDistance, group = ExperimentID), color = "black") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_blank(),  # Remove y-axis title
          legend.title = element_blank(),
          legend.position  = "none") +
    scale_fill_manual(values = color_blind_friendly_colors)
allExperimentsBoxPlot
ggsave(allExperimentsBoxPlot, filename = paste0("~/Keerthana/CellTrajectory/Plots/allexperimentBoxPlot_maxMinNormalisedWithControl.svg"), width = 8, height = 6)
ggsave(allExperimentsBoxPlot, filename = paste0("~/Keerthana/CellTrajectory/Plots/allexperimentBoxPlot_maxMinNormalisedWithControl.png"), width = 8, height = 6)

#Another kind of plot
average_data <- otherComparisonTableNormalised %>%
  group_by(SampleName, ExperimentID) %>%
  summarise(AverageMeasurement = mean(NormalizedAdjustedDistance, na.rm = TRUE)) %>%
  ungroup()
average_data$SampleName <- factor(average_data$SampleName, levels = c("control", "ten", "twenty", "forty"))

overall_average <- average_data %>%
  group_by(SampleName) %>%
  summarise(OverallAverage = mean(AverageMeasurement, na.rm = TRUE)) %>%
  ungroup()

# Plot the averaged data
# Assuming average_data and overall_average data frames are correctly prepared

# Assuming average_data and overall_average are correctly prepared beforehand
a <- ggplot(average_data, aes(x = factor(SampleName), y = AverageMeasurement, color = factor(ExperimentID))) +
  # geom_point(alpha = 0.6) + # Plot points for visual differentiation of data
  geom_line(aes(group = ExperimentID), alpha = 0.4) + # Connect points within each experiment
  scale_color_viridis_d(name = "ExperimentID") +

  # Apply geom_smooth() across all data points to show the overall trend
  geom_smooth(aes(group = 1, color = NULL), # Remove color grouping for the trend line
              method = "lm", se = FALSE, color = "black", size = 1.5, linetype="dashed") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
a
averagePlot <- ggplot() +
  geom_point(data = average_data, aes(x = factor(SampleName), y = AverageMeasurement, color = factor(ExperimentID)), alpha = 0.6) +
  geom_line(data = average_data, aes(x = factor(SampleName), y = AverageMeasurement, group = ExperimentID, color = factor(ExperimentID)), alpha = 0.4) +
  scale_color_viridis_d(name = "ExperimentID") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_blank())+
  # Explicitly specify the group for overall averages
  geom_line(data = overall_average, aes(x = factor(SampleName), y = OverallAverage, group = 1), size = 1.2, color = "black")
# Print the plot
print(averagePlot)
ggsave("~/Keerthana/CellTrajectory/Plots/averagedAllExperiments.svg", a, width = 6, height = 6)
ggsave("~/Keerthana/CellTrajectory/Plots/averagedAllExperimentsLegend.png", a, width = 6, height = 6)

averagedNew <- average_data %>%
  mutate(sample = case_when(
    SampleName == "control" ~ "control",
    TRUE ~ "doublet"
  )) %>%
  select(sample, AverageMeasurement, ExperimentID)

allExperimentsPlotSingletDoublet <- ggplot(averagedNew, aes(x = sample, y = AverageMeasurement)) +
  geom_boxplot(aes(x = sample, y = AverageMeasurement)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_blank(),  # Remove y-axis title
        legend.title = element_blank(),
        legend.position  = "none")
allExperimentsPlotSingletDoublet
ggsave("~/Keerthana/CellTrajectory/Plots/BoxPlot.svg", plot = allExperimentsPlotSingletDoublet, width = 4, height = 6)
#Plot matched lineages for one sample for example
experimentID <- "FM04"
#plot trajectory for lineage 1 for ten, twenty and thirty to compare
trajectoriesPlot <- list()
samples <- c("control", "ten", "twenty", "forty")
inputTrajectoryFolder <- "~/Keerthana/CellTrajectory/data/curvesPlot/"
inputTrajectoryFolder <- "~/Keerthana/singletCode/CellTrajectory/data/curvesPlot/"

for (sample in samples) {
  filepath <- paste0(inputTrajectoryFolder, experimentID, "_", sample, ".csv")
  df <- read_csv(filepath) %>%
    filter(Lineage == 1) %>%
    mutate(Sample = sample)  # Add a column to identify the sample
  trajectoriesPlot[[sample]] <- df
}

# Combine all samples into a single dataframe
combinedTrajectories <- bind_rows(trajectoriesPlot)
combinedTrajectories$Sample <- factor(combinedTrajectories$Sample, 
                                      levels = c("control", "ten", "twenty", "forty"))

combinedTrajectoriesTemp <- combinedTrajectories %>% filter(Sample %in% c("control", "forty"))


# Plotting
plotCurves <- ggplot(combinedTrajectories, aes(x = PC_1, y = PC_2, group = interaction(Sample, Lineage), color = Sample)) +
  geom_path(size = 1.2) + # Draw lines
  scale_color_manual(values = color_blind_friendly_colors) + # Using the color-blind friendly palette
  theme_void() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")


# Display the plot
print(plotCurves)
ggsave(plotCurves, file = "~/Keerthana/CellTrajectory/Plots/FM04_comparison_noLegend.png", width = 4, height = 6)
svglite(filename = "~/Keerthana/CellTrajectory/Plots/FM04_comparison_noLegend.svg", width = 4, height = 6)
plot(plotCurves)
dev.off()
# Filter control and ten samples
control_points <- combinedTrajectoriesTemp %>%
  filter(Sample == "control") %>%
  mutate(match_id = Order)  # Create an identifier to match points

ten_points <- combinedTrajectoriesTemp %>%
  filter(Sample == "ten") %>%
  mutate(match_id = row_number())  # Ensure this matches the control points' order

forty_points <- combinedTrajectoriesTemp %>%
  filter(Sample == "forty") %>%
  mutate(match_id = Order)  # Ensure this matches the control points' order

# Assuming both dataframes are in the same order and have the same number of points
matched_points <- merge(control_points, forty_points, by = "match_id", suffixes = c("_control", "_forty"))

otherColours <- c("hotpink3", "lightblue3")

plotCurves <- ggplot() +
  geom_path(data = combinedTrajectoriesTemp, aes(x = PC_1, y = PC_2, group = interaction(Sample, Lineage), color = Sample)) +
  geom_segment(data = matched_points, 
               aes(x = PC_1_control, y = PC_2_control, xend = PC_1_forty, yend = PC_2_forty),
               color = "black") + # Draw lines between matched points
  geom_point(data = combinedTrajectoriesTemp, aes(x = PC_1, y = PC_2, color = Sample)) +
  scale_color_manual(values = otherColours) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom")

# Display the plot
print(plotCurves)
ggsave(plotCurves, file = "~/Keerthana/CellTrajectory/Plots/FM04_comparison_Illustration_Legend.png", width = 4, height = 6)
svglite(filename = "~/Keerthana/CellTrajectory/Plots/FM04_comparison_Illustration.svg", width = 4, height = 6)
plot(plotCurves)
dev.off()

#Boxplot for just one experiment - FM04
otherComparisonTableNormalised_FM04 <- otherComparisonTableNormalised %>%
  filter(ExperimentID == "FM04")

median_data <- otherComparisonTableNormalised_FM04 %>%
  group_by(ExperimentID, SampleName) %>%
  summarize(MedianNormalizedAdjustedDistance = median(NormalizedAdjustedDistance, na.rm = TRUE)) %>%
  ungroup()

FM04Plot <- ggplot(otherComparisonTableNormalised_FM04, aes(x = SampleName, y = NormalizedAdjustedDistance, fill = SampleName)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  facet_wrap(~ExperimentID, scales = "free_x") +
  geom_line(data = median_data, aes(x = SampleName, y = MedianNormalizedAdjustedDistance, group = ExperimentID), color = "black") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_blank(),  # Remove y-axis title
        legend.title = element_blank(),
        legend.position  = "none") +
  scale_fill_manual(values = color_blind_friendly_colors)
FM04Plot
ggsave(FM04Plot, filename = paste0("~/Keerthana/CellTrajectory/Plots/boxPlotFM04.png"), width = 6, height = 4)
ggsave(FM04Plot, filename = paste0("~/Keerthana/CellTrajectory/Plots/boxPlotFM04.svg"), width = 6, height = 4)
