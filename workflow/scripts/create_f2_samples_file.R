#!/usr/bin/env Rscript

# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Script

library(tidyverse)

# Get list of files
files = list.files(snakemake@input[["input_dirs"]],
                   full.names = T,
                   pattern = "sequence.txt.gz$")

# Process and output samples file
data.frame("PATH" = files) %>%
    # Get basenames
    dplyr::mutate(BASENAME = basename(PATH)) %>% 
    # Get pool numbers
    dplyr::mutate(POOL = stringr::str_split(BASENAME, "_", simplify = T) %>% 
                      subset(select = 3),
                  SAMPLE = stringr::str_split(BASENAME, "_", simplify = T) %>% 
                      subset(select = 6) %>% 
                      # remove "lane" prefix
                      stringr::str_remove(., "lane") %>% 
                      # remove first digit, which is the same as POOL
                      sub(".", "", .) %>% 
                      # conver to integer
                      as.integer(.),
                  PAIR = stringr::str_split(BASENAME, "_", simplify = T) %>% 
                      subset(select = 7)) %>% 
    dplyr::arrange(SAMPLE) %>% 
    # pivot wider to put fq paths into same row
    tidyr::pivot_wider(id_cols = c(SAMPLE, POOL),
                        names_from = PAIR,
                        names_prefix = "fq",
                        values_from = PATH) %>% 
    # write to file
    readr::write_tsv(snakemake@output[[1]])