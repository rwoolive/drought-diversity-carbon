#!/bin/bash


#SBATCH --job-name=rocky_dada_ITS
#SBATCH --output=rocky_dada_ITS.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rwoolive@utk.edu



module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

Rscript ./rocky_dada_ITS/01_dada_ITS.R \
