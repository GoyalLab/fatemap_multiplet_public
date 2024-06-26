#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH --job-name dblce_bcmk
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module purge
module load geos/3.8.1
module load R/4.1.1

cd /projects/p31666/zzhang/doublet-bchmk/repo/benchmarking/scDblFinders
Rscript benchmarking_paper.R
