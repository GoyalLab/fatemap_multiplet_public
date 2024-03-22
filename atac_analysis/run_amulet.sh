#!/bin/bash
#SBATCH --account p31666
#SBATCH --partition genhimem
#SBATCH --job-name amulet
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 600G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


date

# Module prep
module purge all
module load R/4.2.3

cd /projects/p31666/zzhang/doublet-bchmk/repo/atac_analysis
Rscript amulet.R