FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="93814b1d4d0538891a4e9dcb8d36e880a14bce334e6700f7f1d22a92abdbd4de"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https:/github.com/snakemake/snakemake-wrappers/raw/0.74.0/bio/minimap2/aligner/environment.yaml
#   prefix: /conda-envs/9a3940a46748c5fd6f65a34cc26708b7
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - minimap2 ==2.17
RUN mkdir -p /conda-envs/9a3940a46748c5fd6f65a34cc26708b7
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.74.0/bio/minimap2/aligner/environment.yaml /conda-envs/9a3940a46748c5fd6f65a34cc26708b7/environment.yaml

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

RUN mamba env create --prefix /conda-envs/9a3940a46748c5fd6f65a34cc26708b7 --file /conda-envs/9a3940a46748c5fd6f65a34cc26708b7/environment.yaml && \
    mamba env create --prefix /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57 --file /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57/environment.yaml && \
    mamba env create --prefix /conda-envs/9608721699f97513ba7f47bd4e3db24b --file /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml && \
    mamba clean --all -y
