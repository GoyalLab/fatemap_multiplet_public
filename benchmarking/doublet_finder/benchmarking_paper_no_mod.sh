#!/bin/bash
#SBATCH -A p31666
#SBATCH -p short
#SBATCH --job-name nm_df_bcmk
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 1:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/nm_hbrd_bcmk%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/nm_hbrd_bcmk%j.%N.err
#SBATCH --verbose

module purge all
module load R/4.2.3

cd /projects/p31666/zzhang/doublet-bchmk/repo/benchmarking/doublet_finder
Rscript benchmarking_paper_no_mod.R