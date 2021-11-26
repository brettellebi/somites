#IN_DIR = "/hps/nobackup/birney/users/ian/somites/dpABs/F2/no_repeat_reads_or_pers_hets"
#IN_FILES = list.files(IN_DIR, full.names = T)

IN_FILES = unlist(snakemake@input)
OUT_FILE = snakemake@output[[1]]
LOG_FILE = snakemake@log[[1]]

# Send output to log

log <- file(LOG_FILE, open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Read in files

names(IN_FILES) = IN_FILES %>% 
  basename(.) %>% 
  stringr::str_remove(".txt")

df = purrr::map(IN_FILES, function(SAMPLE){
  readr::read_tsv(SAMPLE,
                  col_names = c("CHROM", "POS", "CAB_ALLELE", "CAB", "KAGA_ALLELE", "KAGA"),
                  col_types = c("iicici"))
}) %>% 
  dplyr::bind_rows(.id = "SAMPLE")

out = df %>% 
  dplyr::select(-c(CAB_ALLELE, KAGA_ALLELE)) %>% 
  tidyr::pivot_wider(id_cols = c("CHROM", "POS"),
                     names_from = SAMPLE, values_from = c("CAB", "KAGA"),
                     names_glue = "{SAMPLE}_{.value}")

# Create vector for column order
sorted_names = names(IN_FILES) %>% 
  as.numeric() %>% 
  sort()

col_order = purrr::map(sorted_names, function(SAMPLE) {
  c(paste(SAMPLE, "CAB", sep = "_"),
    paste(SAMPLE, "KAGA", sep = "_"))
}) %>% unlist()


# Order columns and write to file
out %>% 
  dplyr::select(CHROM, POS, all_of(col_order)) %>% 
  # write to file
  readr::write_tsv(OUT_FILE, col_names = T)
