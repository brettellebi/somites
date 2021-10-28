#!/bin/bash

# NOTE: raw data locations:
# FTP: /nfs/ftp/private/indigene_ftp/upload/Ali
# Codon: /hps/nobackup/birney/projects/indigene/raw_data/Ali/

####################
# EBI codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
cd /hps/software/users/birney/ian/repos/somites
conda activate snakemake_6.7.0
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s workflow/Snakefile \
  -p

####################
# RStudio Server
####################

# Build container
# Set container path
CONT=/hps/software/users/birney/ian/containers/somites/R_4.1.0.sif

singularity build --remote \
    $CONT \
    workflow/envs/R_4.1.0/R_4.1.0.def

ssh proxy-codon
bsub -M 50000 -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
CONT=/hps/software/users/birney/ian/containers/somites/R_4.1.0.sif
singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/software/users/birney/ian/tmp:/tmp \
                  --bind /hps/software/users/birney/ian/run:/run \
                  $CONT
# Then run rserver, setting path of config file containing library path
rserver --rsession-config-file /hps/software/users/birney/ian/repos/somites/workflow/envs/rstudio_server/rsession.conf

ssh -L 8787:hl-codon-37-04:8787 proxy-codon

####################
# Copying data from FTP to Codon cluster
####################

# e.g.
bsub -o /dev/null -q datamover "cp -r /nfs/ftp/private/indigene_ftp/upload/Ali/Kaga-Cab_F2_Fish201-400_WGS /nfs/research/birney/projects/indigene/raw_data/Ali/"

# First batch `Kaga-Cab_F2_First200WGS`: 186 - 1 (171) = 185
# Second batch `Kaga-Cab_F2_Fish201-400_WGS`: 192