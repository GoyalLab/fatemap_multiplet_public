# Fatemap Barcode Project

This is the official public code repo of the Fatemap Barcode project. 

The code to perform each doublet detection tool on 
our dataset is saved in the `benchmarking` folder. 

The code to indentify ground truth singlets is deposited in the 
`singlet_ground_truth_counting` folder. 

# Fatemap Data Structure

The data is deposited at this link. It has the following stucture. All the `singlet_pair` files indicate which singlets
are used to simulate doublets in which sample. THe `singlets.txt` file contains all the singlets for that dataset. 

```bash
├── FM01
│   ├── FM01_1_singlet_pairs.csv
│   ├── FM01_2_singlet_pairs.csv
│   ├── ...
│   ├── FM01_singlets.txt
├── FM02
├── FM03
├── ...
└── FM08
```