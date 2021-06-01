FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="9733d02eb40010c8727a5b42313bbd6d5a47070b743dbffd2fed24ffb0962d9e"

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

# Conda environment:
#   source: workflow/envs/minimap2_2.19.yaml
#   prefix: /conda-envs/7294d5ab39de245573a8a3536f1949a0
#   name: minimap2_2.19
#   channels:
#     - bioconda
#     - defaults
#   dependencies:
#     - _libgcc_mutex=0.1=main
#     - k8=0.2.5=he513fc3_0
#     - libgcc-ng=9.1.0=hdf63c60_0
#     - libstdcxx-ng=9.1.0=hdf63c60_0
#     - minimap2=2.17=hed695b0_3
#     - zlib=1.2.11=h7b6447c_3
#   prefix: /hps/software/users/birney/ian/miniconda3/envs/minimap2_2.19
RUN mkdir -p /conda-envs/7294d5ab39de245573a8a3536f1949a0
COPY workflow/envs/minimap2_2.19.yaml /conda-envs/7294d5ab39de245573a8a3536f1949a0/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57 --file /conda-envs/bcb2bf3e95291b3b66344f7fa6fd6a57/environment.yaml && \
    mamba env create --prefix /conda-envs/9608721699f97513ba7f47bd4e3db24b --file /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml && \
    mamba env create --prefix /conda-envs/7294d5ab39de245573a8a3536f1949a0 --file /conda-envs/7294d5ab39de245573a8a3536f1949a0/environment.yaml && \
    mamba clean --all -y
