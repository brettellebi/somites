# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

#PHEN = "/hps/nobackup/birney/users/ian/somites/phens/true/intercept.phen"
#COVAR = "/hps/nobackup/birney/users/ian/somites/covars/true/Microscope.covar"
#MICR = "AU"

## True

PHEN = snakemake@input[["phen"]]
COVAR = snakemake@input[["covar"]]
MICR = snakemake@params[["microscope"]]
OUT_PHEN = snakemake@output[["phen"]]
OUT_LIST = snakemake@output[["ids"]]

# Read in files and combine

phen = readr::read_tsv(PHEN,
                       col_names = c("FID", "IID", "PHEN"))

covar = readr::read_tsv(COVAR,
                        col_names = c("FID", "IID", "Microscope"))

df = left_join(phen, covar, by = c("FID", "IID"))

# Filter for target microscope

out = df %>% 
  dplyr::filter(Microscope == MICR) %>% 
  dplyr::select(FID, IID, PHEN)

# Write to file

readr::write_tsv(out, OUT_PHEN, col_names = F)

# Create IDs list

out = out %>% 
  dplyr::select(FID, IID)

# Write to file

readr::write_tsv(out, OUT_LIST, col_names = F)

