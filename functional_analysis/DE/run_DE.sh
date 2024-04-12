#!/bin/bash
#SBATCH --account p31666
#SBATCH --partition genhimem
#SBATCH --job-name DE
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 500G
#SBATCH --time 8:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


module load R/4.2.3
cd /projects/p31666/zzhang/doublet-bchmk/repo/functional_analysis/DE
Rscript DE.R
