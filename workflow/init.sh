#!/bin/bash

cd /hps/software/users/birney/ian/repos/somites
conda activate snakemake
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -n {cluster.n} -M {cluster.memory} -o {cluster.output} -e {cluster.error}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s workflow/Snakefile \
  -p