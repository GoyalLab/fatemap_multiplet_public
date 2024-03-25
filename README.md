# Fatemap Barcode Project

This is the official public code repo of the Fatemap Barcode project. 

The code to perform each doublet detection tool on 
our dataset is saved in the `benchmarking` folder. 


# Fatemap Data Structure

The data is deposited at this link. It has the following stucture. All the `corrected_singlet_pairs` files indicate which singlets
are used to simulate doublets in which sample. The `singlets.txt` file contains all the singlets for that sample in that dataset. 

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
└── Other Datasets
```

# Singlets definition
We defnie singets based on the following criteria:
1. a single barcode identified per cell ID
2. if there were multiple barcodes per cell, one barcode had significantly more UMI count than other barcodes within the same cell
3. multiple barcodes per cell but the same barcode combination was found in other cells in the same sample
4. multiple barcodes per cell but the same barcode combination was found in other cells across samples within the same experimental design

# Webiste

For more information, please visit our [website]([(https://goyallab.github.io/SingletCodeWebsite]((https://goyallab.github.io/SingletCodeWebsite)) for results, figures, and more!