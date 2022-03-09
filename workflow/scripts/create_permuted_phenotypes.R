# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debugging
#IN_FILE = "/hps/software/users/birney/ian/repos/somites/data/20220214_phenotypes.xlsx"
#PERM_SEED = 1
#OUT_FILE = "/hps/nobackup/birney/users/ian/somites/permuted_phenos/20220214/unsegmented_psm_area/1.xlsx"

## True
IN_FILE = snakemake@input[[1]]
PERM_SEED = snakemake@params[["permutation_seed"]] %>%
  as.numeric()
OUT_FILE = snakemake@output[[1]]

# Read in file and wrangle
phenos = readxl::read_xlsx(IN_FILE)

# Set random seed
set.seed(PERM_SEED)

# Randomise all columns other than `fish`
permute_order = sample(nrow(phenos))
out = phenos %>% 
  dplyr::mutate(dplyr::across(-fish,
                              ~.x[order(permute_order)])
                )

# Write to file
writexl::write_xlsx(out, OUT_FILE)
