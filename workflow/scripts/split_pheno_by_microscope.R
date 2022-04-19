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
OUT = snakemake@output[[1]]

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

readr::write_tsv(out, OUT, col_names = F)
