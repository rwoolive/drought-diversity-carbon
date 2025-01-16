#!/bin/sh


#SBATCH --job-name=ITS_taxon
#SBATCH --output=ITS_taxon.log
#SBATCH --error=ITS_taxon.error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rwoolive@utk.edu




module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2



Rscript ./03_addtax_ITS.R \ \
 

