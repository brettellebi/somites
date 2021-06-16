######################
# Libraries
######################

import pandas as pd
import numpy as np
import os
#from snakemake.utils import validate

#container: "docker://mambaorg/micromamba:0.13.1"
#container: "../../containers/global.sif"
#container: "docker://condaforge/mambaforge:4.9.2-7"

######################
# Config file and sample sheets
######################

configfile: "config/config.yaml"

#validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"], comment = '#', dtype = str).set_index(
    ["sample", "unit"], drop=False
)

F2_samples = pd.read_table(config["F2_samples"], comment = '#', dtype = str).set_index(
    ["LANE"], drop=False
)

with open(config["F2_valid_lanes"]) as file:
    VALID_LANES = [line.strip() for line in file]



#units = pd.read_table(config["units"], dtype=str).set_index(
#    ["sample", "unit"], drop=False
#)
#units.index = units.index.set_levels(
#    [i.astype(str) for i in units.index.levels]
#)  # enforce str in index
#validate(units, schema="../schemas/units.schema.yaml")
########################
## Wildcard constraints
#######################
#
#wildcard_constraints:
#    vartype="snvs|indels",
#    sample="|".join(samples.index),
#    unit="|".join(units["unit"]),
#
#######################
## Helper functions
#######################

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna().tolist()

def get_fastq_F2(wildcards):
    """Get fastq files of given sample-unit."""
    return F2_samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna().tolist()

def get_contigs(start = config["contigs"][0], end = config["contigs"][1]):
    """Get list of chromosomes."""
    end = end + 1
    return list(range(start, end))
