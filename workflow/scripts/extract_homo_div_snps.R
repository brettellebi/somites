# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
GENO_FILE = "/hps/nobackup/birney/users/ian/somites/genos/F0_and_F1/snps_with_AD/all.txt"
MIN_AD = 5

## True
GENO_FILE = snakemake@input[["genos"]]
MIN_AD = snakemake@params[["min_ad"]] %>% 
  as.numeric()
OUT_FULL = snakemake@output[["full"]]
OUT_PASS = snakemake@output[["pass"]]

# Read in file

genos = readr::read_tsv(GENO_FILE)

# Tidy up column names
colnames(genos) = colnames(genos) %>% 
  stringr::str_remove(pattern = "\\[.*\\]") %>% 
  stringr::str_remove(pattern = "# ") %>% 
  stringr::str_replace(pattern = ":", "_")
  
# Recode genotypes

rc = genos %>% 
  # remove allele depth columns for Cab and Kaga
  dplyr::select(-c(Cab_AD, Kaga_AD)) %>% 
  # recode alleles
  dplyr::mutate(dplyr::across(c(Cab_GT,
                                 Kaga_GT,
                                 F1_GT),
                              ~dplyr::case_when(. == "1|1" | . == "1/1" ~ "11",
                                                . == "0|0" | . == "0/0" ~ "00",
                                                . == "0|1" | . == "0/1" ~ "01",
                                                . == ".|." | . == "./." ~ "MISS",
                                                TRUE ~ .),
                              .names = "{.col}_RC")) %>% 
  # split F1 allele depth column
  tidyr::separate(col = F1_AD,
                  into = c("F1_AD_REF", "F1_AD_ALT"),
                  sep = ",",
                  convert = T) %>% 
  # reorder columns
  dplyr::select(CHROM, POS, Cab_GT, Kaga_GT, F1_GT, F1_AD_REF, F1_AD_ALT, Cab_GT_RC, Kaga_GT_RC) %>% 

# Filter for 

filt = rc %>% 
  dplyr::filter(FILTER = dplyr::case_when(Cab_GT == "00" & Kaga_GT == "11" & F1_GT == "01" & F1_REF_AD >= MIN_AD & F1_ALT_AD >= MIN_AD ~ "PASS",
                                          Cab_GT == "11" & Kaga_GT == "00" & F1_GT == "01" & F1_REF_AD >= MIN_AD & F1_ALT_AD >= MIN_AD ~ "PASS",
                                          TRUE ~ "FAIL"))


