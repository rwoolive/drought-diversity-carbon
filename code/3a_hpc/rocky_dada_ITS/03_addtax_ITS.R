## ------------------------------------------------------------------------
#This script adds taxonomy to each of the inferred ASVs

library(dada2)
library(tidyverse)


cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

seqtab.nochim <- "rocky_dada_ITS/output/02_nochimera_mergeruns/rocky_dada_ITS/rocky_dada_ITS_seqtab_final.rds" #change to your file structure
output_ITS <- "rocky_dada_ITS/output"
name <- "rocky_dada_ITS"
tax_db <- "rocky_dada_ITS/sh_general_release_dynamic_27.10.2022.fasta" #change to your tax database (Silva, UNITE, MaarjAM*)
method <- "dada"
threshold <- as.integer(50)

dir.create(file.path(output_ITS, "03_taxonomy"), showWarnings = FALSE)
dir.create(file.path(output_ITS, "03_taxonomy", name), showWarnings = FALSE)

output_ITS <- paste0(output_ITS,"/03_taxonomy/",name,"/")

# Assign taxonomy (general)
set.seed(42) #random  generator necessary for reproducibility

seqtab.nochim <- readRDS(seqtab.nochim)



taxid <- assignTaxonomy(seqtab.nochim, tax_db, multithread = TRUE, tryRC = TRUE)

# Write to disk
saveRDS(taxid, paste0(output_ITS, name, "_tax_assignation.rds"))

cat('\n')
cat(paste0('# The obtained taxonomy file can be found in "', paste0(output_ITS, name, "_tax_assignation.rds"), '"\n'))
cat('\n# All done!\n\n')
