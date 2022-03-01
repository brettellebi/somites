# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
IN_FILES = c("/nfs/research/birney/users/ian/somites/association_testing/20220214/all_sites/permutations/unsegmented_psm_area/5000/1.rds",
             "/nfs/research/birney/users/ian/somites/association_testing/20220214/all_sites/permutations/unsegmented_psm_area/5000/2.rds")

## True
IN_FILES = snakemake@input
OUT_FILE = snakemake@output[[1]]

# Read in files and extract minimum p-value

perm_df = purrr::map_dfr(IN_FILES, function(IN_FILE){
  IN = readRDS(IN_FILE)
  OUT = tibble::tibble(MIN_P = IN$results$p_value_REML %>%
                         min(., na.rm = T)
  )
}, .id = "FILENAME")

# Clean data

perm_df_clean = perm_df %>% 
  # extract variables from file names
  tidyr::separate(col = FILENAME, sep = "/",
                  into = c(rep(NA, 8), "SITE_FILTER", NA, "TARGET_PHENO", "COVARIATES", "INVERSE_NORM", "BIN_LENGTH", "PERMUTATION")) %>% 
  # remove .rds extension from PERMUTATION
  dplyr::mutate(PERMUTATION = stringr::str_remove(PERMUTATION, ".rds"))

# Get minimum
perm_df_mins = perm_df_clean %>% 
  dplyr::group_by(SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, BIN_LENGTH) %>% 
  dplyr::summarise(MIN_P = min(MIN_P))

# Write to file
readr::write_csv(perm_df_mins,
                 perm_file)