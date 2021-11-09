######################
# Libraries
######################

import pandas as pd
import numpy as np
import os

######################
# Config file and sample sheets
######################

configfile: "config/config.yaml"

F0_samples = pd.read_table(config["F0_samples"], comment = '#', dtype = str).set_index(
    ["sample", "unit"], drop=False
)

## Note: THIS RULE WAS RUN BEFORE ANYTHING ELSE,
## WHILE HASHING OUT THE TWO `F2_samples = pd....` lines below
## Run again when config["F2_sequence_dirs"] changes,
## e.g. when this pipeline is run on a new HPC.
rule create_f2_samples_file:
    input:
        input_dirs = config["F2_sequence_dirs"]
    output:
        config["F2_samples"]
    log:
        os.path.join(config["working_dir"], "logs/create_f2_samples_file/create_f2_samples_file.log")
    container:
        "docker://brettellebi/somite_f2:latest"
    script:
        "../scripts/create_f2_samples_file.R"


F2_samples = pd.read_table(config["F2_samples"], comment = '#', dtype = str).set_index(
    ["SAMPLE"], drop=False
)

PERM_SEEDS = list(range(1, config["n_permutations"][0] + 1))

#######################
## Helper functions
#######################

#def get_fastq_F0(wildcards):
#    """Get fastq files of given sample-unit."""
#    return F0_samples.loc[(wildcards.F0_sample, wildcards.unit), ["fq1", "fq2"]].dropna().tolist()

#def get_fastq_F2(wildcards):
#    """Get fastq files of given sample-unit."""
#    return F2_samples.loc[(wildcards.F2_sample), ["fq1", "fq2"]].dropna().tolist()

def get_contigs(start = config["contigs"][0], end = config["contigs"][1]):
    """Get list of chromosomes."""
    end = end + 1
    return list(range(start, end))
