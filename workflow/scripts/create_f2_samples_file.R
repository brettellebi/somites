# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Get variables

## Debug
IN_DIRS = c("/nfs/ftp/private/indigene_ftp/upload/Ali/Kaga-Cab_F2_First200WGS",
            "/nfs/ftp/private/indigene_ftp/upload/Ali/Kaga-Cab_F2_Fish201-400_WGS",
            "/nfs/ftp/private/indigene_ftp/upload/Ali/Kaga-Cab_F2_Fish401-595_WGS")

## True
IN_DIRS = snakemake@input[["input_dirs"]]
OUT_FILE = snakemake@output[[1]]

# Script

library(tidyverse)

# Print variables

print("IN_DIRS:")
print(IN_DIRS)

# Get list of files
files = list.files(IN_DIRS,
                   full.names = T,
                   pattern = "sequence.txt.gz$")

# Process and output samples file
out = data.frame("PATH" = files) %>%
    # Get basenames
    dplyr::mutate(BASENAME = basename(PATH)) %>% 
    # Split to get POOL, SAMPLE, and PAIR variables
    tidyr::separate(col = BASENAME,
                    into = c(NA, NA, "POOL", NA, NA, "SAMPLE", "PAIR", NA),
                    sep = "_") %>% 
    # Tidy up SAMPLE variable
    dplyr::mutate(SAMPLE = SAMPLE %>% 
                      # remove "lane" prefix
                      stringr::str_remove(., "lane.") %>% 
                      as.integer(.)) %>% 
    dplyr::arrange(SAMPLE) %>% 
    # pivot wider to put fq paths into same row
    tidyr::pivot_wider(id_cols = c(SAMPLE, POOL),
                        names_from = PAIR,
                        names_prefix = "fq",
                        values_from = PATH)


# Write to file
readr::write_tsv(out, OUT_FILE)
