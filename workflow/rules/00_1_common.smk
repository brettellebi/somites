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
## WHILE HASHING OUT all modules other than this one
## Run again when config["F2_sequence_dirs"] changes,
## e.g. when this pipeline is run on a new HPC.
rule create_F1_samples_file:
    input:
        input_dirs = config["F1_sequence_dirs"]
    output:
        config["F1_samples"]
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_f1_samples_file/create_f1_samples_file.log"
        )
    # need to use conda env here because of a quirk on the codon cluster at EBI:
    # you can't see the FTP directories while inside a container
    conda:
        "../envs/tidyverse_1.3.1.yaml"
    resources:
        mem_mb = 500
    script:
        "../scripts/create_F1_samples_file.R"

## Note: THIS RULE WAS RUN BEFORE ANYTHING ELSE,
## WHILE HASHING OUT all modules other than this one
## Run again when config["F2_sequence_dirs"] changes,
## e.g. when this pipeline is run on a new HPC.
rule create_F2_samples_file:
    input:
        input_dirs = config["F2_sequence_dirs"]
    output:
        config["F2_samples"]
    log:
        os.path.join(config["working_dir"], "logs/create_f2_samples_file/create_f2_samples_file.log")
    # need to use conda env here because of a quirk on the codon cluster at EBI:
    # you can't see the FTP directories while inside a container
    conda:
        "../envs/tidyverse_1.3.1.yaml"
    resources:
        mem_mb = 500
    script:
        "../scripts/create_F2_samples_file.R"

