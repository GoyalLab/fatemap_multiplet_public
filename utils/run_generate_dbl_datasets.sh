#!/bin/bash
#SBATCH --account p31666
#SBATCH --partition genhimem
#SBATCH --job-name DE_data
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 480G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/DE_data.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/DE_data.err
#SBATCH --verbose

module purge
source /home/zzj4347/.bashrc
module load R/4.2.3

cd /projects/p31666/zzhang/doublet-bchmk/repo/utils
Rscript generate_dbl_datasets.R