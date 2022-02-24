######################
# Libraries
######################

import pandas as pd
import numpy as np
import os

######################
# Get F1 and F2 samples
######################

F1_samples = pd.read_table(config["F1_samples"], comment = '#', dtype = str).set_index(
    ["SAMPLE"], drop=False
)

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
