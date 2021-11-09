# Fetch variables

#IN_FILE = "/nfs/research/birney/users/ian/mikk_genome/repeats/medaka_hdrr_repeats.fixed.gff"
#REPEATS_FILE = "/nfs/research/birney/users/ian/somites/repeats/hdrr_repeats.bed"
IN_FILE = snakemake@input[[1]]
OUT_FILE = snakemake@output[[1]]
LOG_FILE = snakemake@log[[1]]

# Send log

log <- file(LOG_FILE, open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Read in data and write to file

read.table(IN_FILE,
           header = F, sep = "\t", skip = 3, comment.char = "", quote = "", as.is = T) %>%
  dplyr::select(chr = V1, start = V4, end = V5) %>% 
  readr::write_delim(OUT_FILE, delim = "\t", col_names = F)

