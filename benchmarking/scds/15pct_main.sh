#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH --job-name 15pct_hbrd
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 120G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module purge all
module load R/4.2.3

cd /projects/p31666/zzhang/doublet-bchmk/repo/benchmarking/scds
Rscript 15pct_main.R
