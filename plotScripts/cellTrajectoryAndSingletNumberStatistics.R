### cell trajectory statistical analyses for Zhang Melzer et al 2024
## created 20240321
## last updated 20240321



data = read.csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/Revisions/plotData/cellTrajectory/finalTableWithPenalty.csv")

data$singletsOnly = 0

wilcox.test(data$singletsOnly, data$NormalizedAdjustedDistance, paired = FALSE) # W = 0, p-value < 2.2e-16




data = read.csv("/Users/mem3579/Library/CloudStorage/GoogleDrive-madelinemelzer22@gmail.com/.shortcut-targets-by-id/1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl/ZhangMelzerEtAl/Revisions/plotData/totalSingletNumbers/newCompiledNumbers.csv")


wilcox.test(data$numPaperSingletsQC, data$numTrueSingletsQC, paired = TRUE) # V = 1165, p-value = 0.002399

num_rows(data)
sum(data$num10xCellsQC)
numSingletsQC = sum(data$numTrueSingletsQC) # 293618
sum(data$numTrueSinglets) # 304740

numBarcodedCellsQC = sum(data$numBarcodedCellsQC) #338948
sum(data)


recovery = numSingletsQC/numBarcodedCellsQC

barcodedDoubDetect = data %>% filter(ExperimentName %in% c("FateMap_FM01", "FateMap_FM02", "FateMap_FM03", "FateMap_FM04", "FateMap_FM05", "FateMap_FM06", "FateMap_FM08", "JainEtAl", "JiangEtAl", "SPLINTR"))

numSingletCodeSinglets = sum(barcodedDoubDetect$numTrueSingletsQC)
numPaperSinglets = sum(barcodedDoubDetect$numPaperSingletsQC)

pctMoreSinglets = (numSingletCodeSinglets-numPaperSinglets)/(numPaperSinglets)  #0.2945131
