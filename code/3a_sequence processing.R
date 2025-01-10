# Primer pools were used to amplify the V1-V9 region of the 16S rRNA gene (23 primers) 
# and fungal ITS1 gene (15 primers) using the xGen 16S Amplicon Panel v2 protocol 
# (https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/protocol/xgen-16s-amplicon-panel-v2-and-xgen-its1-amplicon-panel-protocol.pdf?sfvrsn=acace007_16; Integrated DNA Technologies, Inc.). 
# Amplicons were sequenced using NextSeq at the University of Tennessee Genomics Sequencing 
# Core facility. We used the R-based DADA2 pipeline to filter and group sequences into 
# ASVs, then taxonomically categorized ASVs using UNITE ITS and RDP 16S databases 
# (Callahan et al. 2016). 


# # get sample names
# names <- list.files("raw-data/sequence/fastqs/", pattern = 'R1_001.fastq.gz', )
# metadat <- data.frame(file.names = names,
#                       submitted.by = rep(NA, length(names)),
#                       library = rep(NA, length(names)),
#                       well = rep(NA, length(names)),
#                       utia = rep(NA, length(names)),
#                       id  = rep(NA, length(names)))
# for(i in 1:length(names)){
#   chars <- strsplit(names[i], "_")
#   metadat$submitted.by[i] <- chars[[1]][1]
#   metadat$library[i] <- chars[[1]][2]
#   metadat$well[i] <- chars[[1]][3]
#   metadat$utia[i] <- chars[[1]][4]
#   metadat$sample.id[i] <- chars[[1]][5]
# }
# write.csv(metadat, "raw-data/sequence/metadat.csv")
# # after this, add actual sample IDs and treatments
metadat <- read.csv("raw-data/sequence/metadat.csv")



##### ITS


# (perform on Rocky because it is computationally intensive)

##### make plots
# connect to UT VPN
# in terminal:
# ssh rwoolive@rocky.nimbios.org # the passphrase is just my macbook password
# if necessary, remove old files using:
# rm -r rocky_dada_ITS/*
# on my computer's terminal, copy the folder from my local Desktop to rocky:
# scp -r /Users/rachelwooliver/Desktop/rocky_dada_ITS rwoolive@rocky.nimbios.org:~/
  
# if haven't done already, on rocky terminal, remove the sample in 2H1 because it has no reads
# rm rocky_dada_ITS/raw/Wooliver_ITS_2H1_Beever240927_S268_R1_001.fastq.gz
# rm rocky_dada_ITS/raw/Wooliver_ITS_2H1_Beever240927_S268_R2_001.fastq.gz
# on rocky terminal, start the run: 
# sbatch rocky_dada_ITS/00_makeplots_ITS.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less rocky_dada_ITS.log
# (then hit q to go back to main terminal page)
# output will be in the raw and output subfolders
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_ITS/output/00_qprofiles/ /Users/rachelwooliver/Desktop/rocky_dada_ITS/output/



##### dada
# on rocky terminal, start the run: 
# sbatch rocky_dada_ITS/01_dada_ITS.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less rocky_dada_ITS.log
# (then hit q to go back to main terminal page)
# output will be in the raw and output subfolders
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_ITS/output/01_errors-output/ /Users/rachelwooliver/Desktop/rocky_dada_ITS/output/
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_ITS/raw/filtered_rocky_dada_ITS/ /Users/rachelwooliver/Desktop/rocky_dada_ITS/raw/
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_ITS/output/01_errors-output/rocky_dada_ITS/ /Users/rachelwooliver/Desktop/rocky_dada_ITS/output/01_errors-output/


##### remove chimeras
# on rocky terminal, start the run: 
# sbatch rocky_dada_ITS/02_chiremmer_ITS.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less rocky_dada_ITS.log
# (then hit q to go back to main terminal page)
# output will be in the raw and output subfolders
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_ITS/output/02_nochimera_mergeruns/rocky_dada_ITS/ /Users/rachelwooliver/Desktop/rocky_dada_ITS/output/02_nochimera_mergeruns/



# use control-D to log out of rocky


















##### 16S


# (perform on Rocky because it is computationally intensive)

##### make plots
# in terminal:
# ssh rwoolive@rocky.nimbios.org # the passphrase is just my macbook password
# if necessary, remove old files using:
# rm -r rocky_dada_16S/*
# on my computer's terminal, copy the folder from my local Desktop to rocky:
# scp -r /Users/rachelwooliver/Desktop/rocky_dada_16S rwoolive@rocky.nimbios.org:~/
# on rocky terminal, start the run: 
# sbatch rocky_dada_16S/00_makeplots_16S.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less rocky_dada_16S.log
# (then hit q to go back to main terminal page)
# output will be in the raw and output subfolders
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_16S/output/00_qprofiles/ /Users/rachelwooliver/Desktop/rocky_dada_16S/output/

##### dada
# on rocky terminal, start the run: 
# sbatch rocky_dada_16S/01_dada_16S.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less rocky_dada_16S.log
# (then hit q to go back to main terminal page)
# output will be in the raw and output subfolders
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_16S/output/01_errors-output/ /Users/rachelwooliver/Desktop/rocky_dada_16S/output/
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_16S/raw/filtered_rocky_dada_16S/ /Users/rachelwooliver/Desktop/rocky_dada_16S/raw/
# scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_16S/output/01_errors-output/rocky_dada_16S/ /Users/rachelwooliver/Desktop/rocky_dada_16S/output/01_errors-output/



##### remove chimeras
# on rocky terminal, start the run: 
# sbatch rocky_dada_16S/02_chiremmer_16S.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less rocky_dada_16S.log
# (then hit q to go back to main terminal page)
# output will be in the raw and output subfolders
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
 # scp -r rwoolive@rocky.nimbios.org:~/rocky_dada_16S/output/02_nochimera_mergeruns/rocky_dada_16S/ /Users/rachelwooliver/Desktop/rocky_dada_16S/output/02_nochimera_mergeruns/
   
# use control-D to log out of rocky








