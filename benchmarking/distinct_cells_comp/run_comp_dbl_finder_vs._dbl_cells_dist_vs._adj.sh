#!/bin/bash
#SBATCH --account p31666
#SBATCH --partition genhimem
#SBATCH --job-name dist_adj_comp
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 500G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module purge all
module load R/4.2.3


cd /projects/p31666/zzhang/doublet-bchmk/repo/benchmarking/distinct_cells_comp
Rscript comp.R
