Guide to plotsfolder- Zhang, Melzer et al. 2024
Created by MEM on 20240529

The majority of scripts used to make plots are in the plotScripts folder of the fatemap_multiplet_public github, unless stated otherwise: https://github.com/GoyalLab/fatemap_multiplet_public 
The .csv files used to make these plots are in the plotData folder of the singletCode figshare: https://doi.org/10.6084/m9.figshare.25478680   

####### MAIN #######

Figure 1D: 
  data: plotData/cellType/cellTypeNumbersCombined.csv
  plotScript: none, manually created bar graphs using the data in Adobe Illustrator. 

Figure 1E: 
  data: plotData/totalSingletNumbers/newCompiledNumbers.csv
  plotScript: plotSingletDoubletComparison.R
  
Figure 1F: 
  data: plotData/ ***Charles
  plotScript: ***Charles
  
Figure 1G: 
  data: plotData/totalSingletNumbers/newCompiledNumbers.csv
  plotScript: plotSingletDoubletComparison.R
  
Figure 1H, 1I: 
  data: plotData/ ***Karun
  plotScript: kkScript1.R, kkScript2.R
  
Figure 2A: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 243-283

Figure 2B: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 293-317

Figure 2C: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 863-891
  
Figure 2D: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, AUPRC: 430-455, AUROC: 463-483

Figure 2E: 
  data: plotData/averagedAndSummedDoublets
  plotScript: heterogeneityCloseDistant.R, lines 205-232
  

Figure 3A:







Figure 4:





Figure 5:




Figure 6



####### SUPPLEMENTARY #######

Figure S1, S2, S3:
  data: plotData/cellType/UMAPData/
  plotScript: plotScripts/CellTypePlot.R
  
  
Figure S4B:
  data: 
  plotScript: 
  
Figure S4D:
  data: 
  plotScript: 
  
Figure S5: ***Karun 
  data: 
  plotScript: 

Figure S6: 
  data: 
  plotScript: 

Figure S8A: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, AUPRC: 491-528, AUROC: 491-503, 536-556
  
Figure S8B: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, AUPRC: 592-616
  
Figure S11A:
  data: plotData/averagedAndSummedDoublets
  plotScript: heterogeneityCloseDistant.R, lines 205-221, 235-244

Figure S11B:
  data: plotData/averagedAndSummedDoublets
  plotScript: heterogeneityCloseDistant.R, lines 205-221, 248-257
