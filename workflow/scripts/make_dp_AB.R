# Fetch variables

#DP4_FILE = "/hps/nobackup/birney/users/ian/somites/dp4s/F0/all_sites/Cab.dp4.txt"
#SITES_FILE = "/hps/nobackup/birney/users/ian/somites/data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"
DP4_FILE = snakemake@input[["dp4"]]
SITES_FILE = snakemake@input[["sites_file"]]
OUT_FILE = snakemake@output[[1]]
LOG_FILE = snakemake@log[[1]]

# Send log

log <- file(LOG_FILE, open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

#
