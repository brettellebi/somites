#!/bin/bash

# NOTE: raw data locations:
# FTP: /nfs/ftp/private/indigene_ftp/upload/Ali
# Codon: /hps/nobackup/birney/projects/indigene/raw_data/Ali/

####################
# EBI codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -Is bash
#bsub -q datamover -M 20000 -Is bash # when needing to copy raw data from FTP
cd /hps/software/users/birney/ian/repos/somites
conda activate snakemake_6.12.1
# 
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

# To restart jobs with more memory (e.g. for rule run_gwls)
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  --restart-times 0 \
  -s workflow/Snakefile \
  -p

####################
# Containers
####################

HMMCONT=/hps/nobackup/birney/users/ian/containers/somites/hmmlearn_0.2.7.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $HMMCONT \
    workflow/envs/hmmlearn_0.2.7.def

####################
# RStudio Server
####################

# Build general R container
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -q short -Is bash
cd /hps/software/users/birney/ian/repos/somites
CONT=/hps/nobackup/birney/users/ian/containers/somites/R_4.1.0.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $CONT \
    workflow/envs/R_4.1.0/R_4.1.0.def

# Newer R version
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
cd /hps/software/users/birney/ian/repos/somites
CONT=/hps/nobackup/birney/users/ian/containers/somites/R_4.1.3.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $CONT \
    workflow/envs/R_4.1.3/R_4.1.3.def

# Build R container for PhenotypeSimulator
CONT=/hps/nobackup/birney/users/ian/containers/somites/R_4.1.0_PhenotypeSimulator.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $CONT \
    workflow/envs/R_4.1.0/R_4.1.0_PhenotypeSimulator.def

# Run container
ssh proxy-codon
bsub -M 50000 -Is bash
cd /hps/software/users/birney/ian/repos/somites
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
CONT=/hps/nobackup/birney/users/ian/containers/somites/R_4.1.0.sif
CONT=/hps/nobackup/birney/users/ian/containers/somites/R_4.1.3.sif
singularity shell --bind /hps/nobackup/birney/users/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/nobackup/birney/users/ian/tmp:/tmp \
                  --bind /hps/nobackup/birney/users/ian/run:/run \
                  $CONT
# Then run rserver, setting path of config file containing library path
rstudio-server kill-all
rserver \
    --rsession-config-file /hps/software/users/birney/ian/repos/somites/workflow/envs/R_4.1.3/rsession.conf \
    --server-user brettell


ssh -L 8787:hl-codon-37-04:8787 proxy-codon

####################
# Copying data from FTP to Codon cluster
####################

# e.g.
bsub -o /dev/null -q datamover "cp -r /nfs/ftp/private/indigene_ftp/upload/Ali/Kaga-Cab_F2_Fish201-400_WGS /nfs/research/birney/projects/indigene/raw_data/Ali/"

# First batch `Kaga-Cab_F2_First200WGS`: 186 - 1 (171) = 185
# Second batch `Kaga-Cab_F2_Fish201-400_WGS`: 192