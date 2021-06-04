######################
# Libraries
######################

import pandas as pd
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
#
def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna().tolist()
#    if len(fastqs) == 2:
#        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
#    return {"r1": fastqs.fq1}
#

def get_contigs(start = config["contigs"][0], end = config["contigs"][1]):
    """Get list of chromosomes."""
    end = end + 1
    return list(range(start, end))
#def is_single_end(sample, unit):
#    """Return True if sample-unit is single end."""
#    return pd.isnull(units.loc[(sample, unit), "fq2"])
#
#def get_trimmed_reads(wildcards):
#    """Get trimmed reads of given sample-unit."""
#    if not is_single_end(**wildcards):
#        # paired-end sample
#        return expand(
#            os.path.join(config["working_dir"], "fastqs/trimmed/{sample}-{unit}.{group}.fastq.gz"),
#            group=[1, 2],
#            **wildcards
#        )
#    # single end sample
#    return os.path.join(config["working_dir"], "fastqs/trimmed/{sample}-{unit}.{group}.fastq.gz").format(**wildcards)