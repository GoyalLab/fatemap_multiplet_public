# singletCode

This is the official public code repo of analysis scripts for our paper titled: **Synthetic DNA barcodes identify singlets in scRNA-seq datasets and evaluate doublet algorithms**.

This folder contains various analysis scripts used throughout the paper. Some key folders are:
- `benchmarking`: scripts to benchmark the tools on barcoded datasets
- `functional_analysis`: scripts to assess doublets' impact on downstream analysis
- `atac_analysis`: scripts to assess Amulet performance on barcode singlets
- `plotScripts`: scripts used for plotting. INSTRUCTIONS TO REPRODUCE ALL PLOTS FROM PLOTDATA (figshare) ARE IN guideToPlotScripts.md
- `heterogeneity`: scripts to assess heterogeneity across samples.
- `otherBarcodingMethods`: scripts to adapt data from other barcoding technologies to have singlets identified with singletCode.
- `cellTypeAnnotation`: scripts to annotate cell type.
- `cellCount`: scripts to process data according to the original publications for the purpose of comparing cell numbers at different stages of processing, including singlet recovery.
  
# Data Structure

The datasets all have the following structure. All the `corrected_singlet_pairs` files indicate which singlets
are used to simulate doublets in which sample. The `singlets_all.txt` file contains all the singlets for that sample in that dataset. 

```bash
├── FM01
│   ├── 10X   
│   │   ├── sample1
│   │   │   ├── corrected_singlet_pairs.csv
│   │   │   ├── singlets_all.txt
│   │   │   ├── matrix.mtx.gz
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   ├── ...
│   │   ├── sample4
├── FM02
├── FM03
├── ...
├── ClonMapper
├── SPLINTR
├── TREX
└── Other Datasets
```

# Singlets definition
We defnie singets based on the following criteria:
1. a single barcode identified per cell ID
2. if there were multiple barcodes per cell, one barcode had significantly more UMI count than other barcodes within the same cell
3. multiple barcodes per cell but the same barcode combination was found in other cells in the same sample
4. multiple barcodes per cell but the same barcode combination was found in other cells across samples within the same experimental design

# Webiste

For more information, please visit our [website](https://goyallab.github.io/SingletCodeWebsite) for results, figures, and more!
