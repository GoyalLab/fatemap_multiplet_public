# Fatemap Barcode Project

This is the official public code repo of the Fatemap Barcode project. 

The code to perform each doublet detection tool on 
our dataset is saved in the `benchmarking` folder. 

The code to indentify ground truth singlets is deposited in the 
`singlet_ground_truth_counting` folder. 

# Fatemap Data Structure

The data is deposited at this link. It has the following stucture. 

```bash
├── FM01
│   ├── 10X
│   │   ├── sample1
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   ├── sample2
│   │   ├──  ...
│   │   ├── sampleX
│   ├── fatemapID
│   │   ├── FM01_good_data.tsv
│   │   ├── FM01_multiplets.txt
│   │   ├── FM01_singlets.txt
│   │   ├── stepFourStarcodeShavedReads50.txt
├── FM02
├── FM03
├── ...
└── FM08
```