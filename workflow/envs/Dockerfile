FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7045314d10d02b64d481ef3139e5807b4d80c3b21a396e5bcfca7a0392f09e3e"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https:/github.com/snakemake/snakemake-wrappers/raw/0.74.0/bio/reference/ensembl-sequence/environment.yaml
#   prefix: /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57
#   channels:
#     - conda-forge
#   dependencies:
#     - curl
RUN mkdir -p /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.74.0/bio/reference/ensembl-sequence/environment.yaml /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57/environment.yaml

# Conda environment:
#   source: https:/github.com/snakemake/snakemake-wrappers/raw/0.74.0/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/9608721699f97513ba7f47bd4e3db24b
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - samtools ==1.10
RUN mkdir -p /conda-envs/9608721699f97513ba7f47bd4e3db24b
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.74.0/bio/samtools/faidx/environment.yaml /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57 --file /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57/environment.yaml && \
    mamba env create --prefix /conda-envs/9608721699f97513ba7f47bd4e3db24b --file /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml && \
    mamba clean --all -y
