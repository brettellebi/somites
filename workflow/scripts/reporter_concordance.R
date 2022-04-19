# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#BIN_LENGTH = as.numeric("5000")
#GENO = "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn_true/None/5000/1.csv"
#COV = "1"
#REPORTER = "16:28706898-28708417"
#PHENO = here::here("data/20220321_phenotypes.xlsx")
#LOW_COV_SAMPLES = c(26, 89, 166, 178, 189, 227, 239, 470, 472, 473, 490, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511)
#MAX_READS = "None"

## True
GENO = snakemake@input[["geno"]]
PHENO = snakemake@input[["pheno"]]
COV = snakemake@params[["cov"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
REPORTER = snakemake@params[["reporter_loc"]]
LOW_COV_SAMPLES = snakemake@params[["low_cov_samples"]]
OUT = snakemake@output[[1]]


# Read in files

## Genos

genos = readr::read_csv(GENO,
                        col_types = "iiiiiidi") %>% 
  # add key variables
  dplyr::mutate(BIN_START = (BIN * BIN_LENGTH) + 1,
                BIN_END = ((BIN + 1) * BIN_LENGTH)) 

## Phenos

phenos = readxl::read_xlsx(PHENO)

## Reporter
rep = c(stringr::str_split(REPORTER, ":", simplify = T)[1],
        stringr::str_split(REPORTER, ":", simplify = T)[2] %>% 
          str_split(., "-", simplify = T) %>% .[,1],
        stringr::str_split(REPORTER, ":", simplify = T)[2] %>% 
          str_split(., "-", simplify = T) %>% .[,2]) %>% 
  as.numeric()
names(rep) = c("CHROM", "START", "END")

# Get concordance

out = genos %>% 
  dplyr::filter(CHROM == rep["CHROM"] & BIN_START < rep["START"] & BIN_END > rep["END"]) %>% 
  # bind with pheno
  dplyr::left_join(phenos %>% 
                     dplyr::select(SAMPLE, REPORTER_PHENO = reporter_pheno),
                   by = "SAMPLE") %>% 
  # recode reporter
  dplyr::mutate(REPORTER_PHENO = dplyr::recode(REPORTER_PHENO,
                                               `-1` = 0,
                                               `0` = 1,
                                               `1` = 2)) %>% 
  # remove low-coverage samples
  dplyr::filter(!SAMPLE %in% LOW_COV_SAMPLES) %>% 
  # remove samples with NA for reporter pheno
  dplyr::filter(!is.na(REPORTER_PHENO))

conc = tibble(COV = COV,
              CONC = sum(out$STATE == out$REPORTER_PHENO, na.rm = T) / nrow(out))

# Write to file

readr::write_csv(conc, OUT)
