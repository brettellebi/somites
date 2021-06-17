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

samples = pd.read_table(config["samples"], comment = '#', dtype = str).set_index(
    ["sample", "unit"], drop=False
)

## Note: THIS RULE WAS RUN BEFORE ANYTHING ELSE,
## WHILE HASHING OUT THE TWO `F2_samples = pd....` lines below
## Run again when the location of the raw sequencing data changes,
## e.g. when this pipeline is run on a new HPC.
rule create_f2_samples_file:
    input:
        input_dir = os.path.join(config["raw_data_dir"], "Kaga-Cab_F2_First200WGS")
    output:
        config["F2_samples"]
    container:
        "docker://brettellebi/somite_f2:latest"
    script:
        "../scripts/create_f2_samples_file.R"


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
