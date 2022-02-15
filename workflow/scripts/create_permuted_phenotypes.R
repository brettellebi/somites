# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debugging
IN_FILE = "/hps/software/users/birney/ian/repos/somites/data/20220214_phenotypes.xlsx"
TARGET_PHENO = "unsegmented_psm_area"
PERM_SEED = 1
OUT_FILE = "/nfs/research/birney/users/ian/somites/permuted_phenos/unsegmented_psm_area/1.xlsx"

## True
IN_FILE = snakemake@input[[1]]
TARGET_PHENO = snakemake@params[["target_phenotype"]]
PERM_SEED = snakemake@params[["permutation_seed"]] %>%
  as.numeric()
OUT_FILE = snakemake@output[[1]]

# Set seed
set.seed(PERM_SEED)

## Read in file and wrangle
phenos = readxl::read_xlsx(IN_FILE) %>%
    # select key columns
    dplyr::select(fish, all_of(TARGET_PHENO)) %>%
    # ensure that the phenotype column is numeric
    dplyr::mutate(dplyr::across(all_of(TARGET_PHENO),
                                ~ as.numeric(.x))) %>%
    # randomise phenotype
    dplyr::mutate(dplyr::across(all_of(TARGET_PHENO),
                                ~ sample(.x)))

# Write to file
writexl::write_xlsx(phenos, OUT_FILE)
