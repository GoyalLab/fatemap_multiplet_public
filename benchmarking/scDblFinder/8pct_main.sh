#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH --job-name 8pct_DBLCELL
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module purge all
module load R/4.2.3

cd /projects/p31666/zzhang/doublet-bchmk/repo/benchmarking/doublet_cells
Rscript 8pct_main.R
