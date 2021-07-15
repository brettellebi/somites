#!/bin/bash

# NOTE: raw data locations:
# FTP: /nfs/ftp/private/indigene_ftp/upload/Ali
# Codon: /hps/nobackup/birney/projects/indigene/raw_data/Ali/

####################
# EBI codon cluster
####################

#module load singularity-3.5.3-gcc-9.3.0-o6v53jz
#module load singularity-3.6.4-gcc-9.3.0-yvkwp5n
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -Is bash
cd /hps/software/users/birney/ian/repos/somites
conda activate snakemake_6.4.1
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -outdir {cluster.outdir} -o {cluster.outfile} -e {cluster.error}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s workflow/Snakefile \
  -p

####################
# RStudio Server
####################

ssh proxy-codon
bsub -M 20000 -Is """
singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/software/users/birney/ian/tmp:/tmp \
                  --bind /hps/software/users/birney/ian/run:/run \
                  docker://brettellebi/somite_f2:latest
"""
# Then run rserver, setting path of config file containing library path
rserver --rsession-config-file /hps/software/users/birney/ian/repos/somites/workflow/envs/rstudio_server/rsession.conf