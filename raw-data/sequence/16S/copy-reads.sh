#!/bin/bash

#SBATCH --mail-user=tfreem10@vols.utk.edu
#SBATCH --mail-type=END,FAIL


#######################################
# Constants
#######################################
# Array of reads files in fastq format
readonly reads=(/lustre/isaac/proj/UTK0204/UTKGenomicsCoreSequencingData/240927_VL00838_4_AAG5YJNM5/UTK0204_Beever_Amplicon_240927/fastq/Wooliver*.fastq.gz)
# Directory to copy the reads into
readonly destination_dir_reads="fastqs"

#######################################
# Main program.
#######################################
mkdir -p "${destination_dir_reads}"

cp --verbose "${reads[@]}" "${destination_dir_reads}"
