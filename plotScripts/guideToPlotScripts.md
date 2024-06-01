Guide to plots folder- Zhang, Melzer et al. 2024
Created by MEM on 20240531

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
  data: plotData/cellType/singlets_by_category_pct.csv
  plotScript: visualization/barplot_celltype.R, lines 61-66
  
Figure 1G: 
  data: plotData/totalSingletNumbers/newCompiledNumbers.csv
  plotScript: plotSingletDoubletComparison.R
  
Figure 1H, 1I: 
  data: plotData/Gini_Plot_data.tsv
  plotScript: giniAnalysisKK/plotGini.R
  
Figure 2A: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 242-282

Figure 2B: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 293-317

Figure 2C: 
  data: /plotData/benchmarking/TNR_plotted_filtered.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 877-904
  
Figure 2D: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, AUPRC: 430-455, AUROC: 463-483

Figure 2E: 
  data: plotData/averagedAndSummedDoublets.csv
  plotScript: heterogeneityCloseDistant.R, lines 205-232
  
Figure 3A, phenotypic volume:
  data: plotData/acrossSampleHeterogeneity/rankedPhenotypicVolumeForPlot.csv
  plotScript: phenotypicVolumePlot.R, lines 198-222
  
Figure 3A, E-distance:
  data: plotData/acrossSampleHeterogeneity/rankedEDistanceForPlot.csv
  plotScript: EDistancePlots.R, lines 238-262

Figure 3A, Shannon diversity:
  data: plotData/acrossSampleHeterogeneity/rankedShannonClusterDiversityForPlot.csv
  plotScript: shannonEntropyPlots.R, lines 176-204

Figure 3A, fU-fL:
  data: plotData/acrossSampleHeterogeneity/rankedDEHeterogeneityForPlot.csv
  plotScript: DEHeterogeneityPlots.R, lines 244-273

Figure 3A, pUp + pDown:
  data: plotData/acrossSampleHeterogeneity/rankedDEHeterogeneityForPlot.csv
  plotScript: DEHeterogeneityPlots.R, lines 191-220
  
Figure 3C:
  data: plotData/heterogenityAnalysis/all.tsv
  plotScript: heterogeneityCloseDistant.R, lines 88-107
  
Figure 3D:
  data: plotData/heterogenityAnalysis/all.tsv
  plotScript: heterogeneityCloseDistant.R, lines 110-120

Figure 4A:
  data: plotData/similarityScore/mean_similarityMatrix.csv
  plotScript: overlappingSingletsIdentified.R, lines 321-346
  
Figure 4B:
  data: plotData/similarityScore/allPoints_similarityScores.csv
  plotScript: heterogeneityCloseDistant.R, lines 142-176
  
Figure 4C:
  data: plotData/ensembleMethods/chordHybridEnsemblePlotted.csv
  plotScript: ensembleMethodsPlot.R, lines 314-346
  
Figure 4D:
  data: plotData/ensembleMethods/chordHybridEnsemblePlotted.csv
  plotScript: ensembleMethodsPlot.R, lines 352-376
  
Figure 4E: 
  data: plotData/benchmarking/barcodedAllForPlot_ss3.csv
  plotScript: smartseq3Plots.R, lines 87-108
  
Figure 4F: 
  data: plotData/benchmarking/barcodedAllForPlot_ss3.csv
  plotScript: smartseq3Plots.R, lines 115-134
  
Figure 4H: 
  data: plotData/AmuletAnalysis/combined_umi1_res.csv
  plotScript: heterogeneityCloseDistant.R, lines 181-201

Figure 5B:
  data: plotData/differentialExpression/de.csv
  plotScript: differentialExpressionPlot.R
  
Figure 5C:
  data: plotData/clusteringStability/resultsclusteringStability_Louvain_allDatasets.csv
  plotScript: clusteringStabilityPlot.R, lines 76-102, 114-127

Figure 5D:
  data: plotData/cellCommunication/cc_formatted.csv
  plotScript: cellChatPlot.R
  
Figure 5E:
  plotScript: functional_analysis/CellTrajectory/SlingShot2.R
  
Figure 5G:
  data: plotData/cellTrajectory/finalTableWithPenalty.csv
  plotScript: functional_analysis/CellTrajectory/graphDistance.R, line plot: lines 260-299, 319-329, 348-360, bar plot: 364-379

Figure 6B: 
  data: plotData/classifier/allClassifierData_formatted.csv
  plotScript: classifierPlots.R, lines 249-273
  
Figure 6C: 
  data: plotData/classifier/allClassifierData_formatted.csv
  plotScript: classifierPlots.R, lines 300-334
  
Figure 6D: 
  data: plotData/classifier/allClassifierData_formatted.csv
  plotScript: classifierPlots.R, lines 249-250, 277-296

Figure 6E: 
  data: plotData/classifier/allClassifierData_formatted.csv
  plotScript: classifierPlots.R, lines 300-301, 339-358
  
Figure 6F: 
  data: plotData/classifier/allClassifierVariableData.csv, plotData/benchmarking/barcodedAllForPlot.csv
  plotScript: classifierPlots.R, lines 394-470 (AUPRC), lines 477-494 (AUROC)
  
Figure 6H: 
  data: plotData/classifier/alls1s2ForPlot_formatted.csv
  plotScript: classifierPlots.R, lines 540, 577-606
  
Figure 6I: 
  data: plotData/classifier/alls1s2ForPlot_formatted.csv
  plotScript: classifierPlots.R, lines 540, 612-627

####### SUPPLEMENTARY #######

Figure S1, S2, S3:
  data: plotData/cellType/UMAPData/
  plotScript: plotScripts/CellTypePlot.R
  
Figure S4B:
  data: plotData/neighborAnalysis/neighborAnalysisPlotData.csv
  plotScript: neighborAnalysis/neighborAnalysis.R, lines 1-166, 258-399
  
Figure S4D:
  data: plotData/neighborAnalysis/neighborAnalysisPlotData.csv
  plotScript: neighborAnalysis/neighborAnalysis.R, lines 220-244
  
Figure S5:
  data: 10X matrices
  plotScript: giniAnalysisKK/extractGini.R

Figure S6: 
  data: plotData/totalSingletNumbers/dataCounts/ and /dataUMAP/
  plotScript: plot10xAll.R
  
Figure S7A:
  data: plotData/benchmarking/calls_barcoded_averageDoublets.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, lines 1291-1323

Figure S7B:
  data: plotData/heatmap/
  plotScript: vizualization/heatmap.py
  
Figure S7C:
  data: plotData/heatmap/
  plotScript: vizualization/heatmap.py

Figure S8A: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, AUPRC: 491-528, AUROC: 491-503, 536-556
  
Figure S8B: 
  data: plotData/benchmarking/barcodedNonBarcoded_AUPRC_AUROC_TNR.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 596-620
  
Figure S9A, phenotypic volume:
  data: plotData/acrossSampleHeterogeneity/rankedPhenotypicVolumeForPlot_AUPRC.csv
  plotScript: phenotypicVolumePlot.R, lines 290-314
  
Figure S9A, E-distance:
  data: plotData/acrossSampleHeterogeneity/rankedEDistanceForPlot_AUROC.csv
  plotScript: EDistancePlots.R, lines 330-354

Figure S9A, Shannon diversity:
  data: plotData/acrossSampleHeterogeneity/rankedShannonClusterDiversityForPlot_AUROC.csv
  plotScript: shannonEntropyPlots.R, lines 268-292

Figure S9A, fU-fL:
  data: plotData/acrossSampleHeterogeneity/rankedDEHeterogeneityForPlot_AUROC.csv
  plotScript: DEHeterogeneityPlots.R, lines 401-428

Figure S9A, pUp + pDown:
  data: plotData/acrossSampleHeterogeneity/rankedDEHeterogeneityForPlot_AUROC.csv
  plotScript: DEHeterogeneityPlots.R, lines 342-368
  
Figure S10A: 
  data: plotData/averagedAndSummedDoublets/averagedAndSummedDoublets.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 1002-1029
  
Figure S10B: 
  data: plotData/averagedAndSummedDoublets/averagedAndSummedDoublets.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 1039-1063

Figure S10C: 
  data: plotData/averagedAndSummedDoublets/TNR_plotted_formatted_sum.csv
  plotScript: benchmarkingResultPlots.R, lines 20-43, 1128-1155

Figure S11A:
  data: plotData/averagedAndSummedDoublets
  plotScript: heterogeneityCloseDistant.R, lines 205-221, 235-244

Figure S11B:
  data: plotData/averagedAndSummedDoublets
  plotScript: heterogeneityCloseDistant.R, lines 205-221, 248-257
  
Figure S12:
  data: plotData/similarityScore/ (all _similarityMatrix.csv)
  plotScript: overlappingSingletsIdentified.R, lines 361-406
  
Figure S13A: 
  data: plotData/benchmarking/barcodedAllForPlot_ss3.csv
  plotScript: smartseq3Plots.R, lines 87-108
  
Figure S13B: 
  data: plotData/benchmarking/barcodedAllForPlot_ss3.csv
  plotScript: smartseq3Plots.R, lines 115-134
  
Figure S14A:
  data: plotData/cellCommunication/cpdb_formatted.csv
  plotScript: cellChatPlot.R
  
Figure S14B:
  data: plotData/clusteringStability/resultsclusteringStability_Louvain_allDatasets.csv
  plotScript: clusteringStabilityPlot.R, lines 76-102

Figure S15:
  plotScript: functional_analysis/CellTrajectory/SlingShot2.R
  
Figure S16:
  plotScript: classifier/testingClassifier.ipynb
  
Figure S17A:
  data: plotData/classifier/summary_s1s2_bothClassifiers_andSelf.csv
  plotScript: classifierPlots.R, lines 667-695

Figure S17B:
  data: plotData/classifier/summary_s1s2_bothClassifiers_andSelf.csv
  plotScript: classifierPlots.R, lines 667-695
  
Figure S17C:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 602-614, 749-775 
  
Figure S17D:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 619-628, 749-775
  
Figure S17E:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 290-307
  
Figure S17F:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 312-321

Figure S17G:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 329-341
  
Figure S17H:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 346-355
  
Figure S17I:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 680-692, 819-831
  
Figure S17J:
  data: plotData/classifier/crossExperiment/
  plotScript: classifierCrossExperimentPlots.R, lines 697-706, 836-845