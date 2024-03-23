# neighborAnalysis for doublet project
Created 20230707 by Madeline E Melzer
Last Modified 20230710 by Madeline E Melzer

This repo contains scripts needed to produce all graphs contained
in figures for the neighbor analyis portion of the paper: Zhang et al., 2023

For any questions, please contact MEM (madeline.melzer at northwestern.edu) and the corresponding
authors for this paper, ZZ ( at northwestern.edu) and YG (yogesh.goyal at northwestern.edu)

In order to reproduce all graphs and images in this analysis:

0. Download and install dependencies for this code:
- R v4.2.2
- R packages tidyverse 2.0.0, Seurat 4.3.0, ggplot2 3.4.2, reshape2 1.4.4, svglite 2.1.1, installr 0.23.4, patchwork 1.1.2, RANN 2.6.1, DescTools 0.99.49, tools 4.2.2, and their associated dependencies

1. Download all extractedData (6.3GB of intermediate-processed sequencing, network, and reprogramming experiment data) files in the structure provided for this repository from OneDrive: 
  https://nuwildcat-my.sharepoint.com/:f:/g/personal/mem3579_ads_northwestern_edu/ElcvvhshalVGtELsug-QuAwBZT0qo-l85z-PoJuQWx9mUQ?e=uZHwEw 

- your top-level project directory should have subdirectories: data and plots.

2. Use quality control thresholds for each dataset indicated in the supplemental information of Zhang et al. 2023 
