# gini analysis for clustering stability to see imbalance of doublet clustering after 


library(tidyverse)
library(dplyr)
library(Seurat)
library(glue)
library(DescTools)
library(ggplot2)

set.seed(23)

resultsDirectory = "R:/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/ZhangMelzerEtAl/data/clusteringStability/results/"

datasetFiles <- list.files(resultsDirectory, pattern = "\\.csv$", full.names = FALSE, recursive = FALSE)

##########################################################################################################################
# gini
##########################################################################################################################
finalResults <- tibble(dataset = character(), dbl_act = numeric(), resolution = numeric(), gini = numeric(), type = character())

for (datasetFile in datasetFiles) {
  
  data <- read.csv(glue("{resultsDirectory}{datasetFile}"))
  
  minMax <- data %>%
    group_by(dbl_act, resolution) %>%
    summarise(min_pct = min(pct_doublets, na.rm = TRUE),
              max_pct = max(pct_doublets, na.rm = TRUE), .groups = 'drop')
  
  # Join this back to the data to have min and max available for each row
  data <- data %>%
    left_join(minMax, by = c("dbl_act", "resolution"))
  
  results <- data %>%
    filter(pct_doublets != 0) %>%
    group_by(dbl_act, resolution) %>%
    summarise(gini = Gini(pct_doublets), .groups = 'drop') %>%
    mutate(type = "actual")
  
  # Add the dataset name to the results
  results$dataset <- first(data$dataset)
  
  # Bind the results for this file to the finalResults dataframe
  finalResults <- bind_rows(finalResults, results)
  
  # Simulate a uniform distribution for each group, calculate the Gini coefficient, and label as "Uniform"
  uniformResults <- data %>%
    filter(pct_doublets != 0) %>%
    group_by(dbl_act, resolution) %>%
    summarise(gini = Gini(runif(n(), min = first(min_pct), max = first(max_pct))), .groups = 'drop') %>%
    mutate(type = "uniform")
  
  # Add the dataset name to the uniformResults
  uniformResults$dataset <- first(data$dataset)
  
  # Bind the uniform distribution results to the finalResults dataframe
  finalResults <- bind_rows(finalResults, uniformResults)
  
  print(results)
  print(uniformResults)
  
}

write.csv(finalResults, glue("{resultsDirectory}clusteringStabilityGini.csv"))



##########################################################################################################################
# plotting
##########################################################################################################################

plotResults = finalResults %>% filter(dbl_act == 40, resolution == 0.2)

plot = ggplot(plotResults, aes(x = gini, fill = type)) +
  geom_histogram(binwidth = 0.025) +
  #labs(x = "Gini Coefficient", y = "Frequency", title = "Distribution of Gini Coefficients: Actual vs Uniform") +
  scale_fill_manual(values = c("actual" = "#8cc53f", "uniform" = "grey85")) +
  facet_wrap(~ type, nrow = 2, scales = "free_y") +
  theme_classic()
plot

  
