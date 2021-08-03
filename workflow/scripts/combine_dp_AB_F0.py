#!/usr/bin/env python3

# Send stdout and stderr to log file

import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

# Import libraries

import pandas as pd

# Read files into list

dfs = list()
for file in snakemake.input:
    df = pd.read_csv(file, sep = "\t", header = None)
    dfs.append(df)

# Bind together

full_df = pd.concat(dfs)

# Write to file

full_df.to_csv(snakemake.output[0], sep = "\t", header = False, index = False)
