# Containerized workflow for calling medaka F2 haplotypes

**NOTE**: Caution re: this issue: https://github.com/snakemake/snakemake/issues/304

## Steps

1. Change directory to your preferred working directory on the cluster, e.g.:
```bash
cd path/to/working/directory
```

2. Clone this repository:
```bash
git clone https://github.com/brettellebi/somites.git
```

3. Install miniconda3:
    - Download appropriate by installer by copying the link from here: https://docs.conda.io/en/latest/miniconda.html, then `wget` it, e.g.:
```bash
wget {link}
```
    - Then follow instructions here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

4. Install Snakemake
