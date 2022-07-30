# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#PED = "/hps/nobackup/birney/users/ian/somites/peds/F2/hdrr/None/5000/1/unsegmented_psm_area.ped"
#PHENO_FILE = here::here("config/DF_648_periodNORMALISED.csv")
#PHENOTYPE = "intercept"

## True
PED = snakemake@input[["ped"]]
PHENO_FILE = snakemake@input[["phenos"]]
PHENOTYPE = snakemake@params[["phenotype"]]
OUT = snakemake@output[[1]]

# Read in files

pedp = readr::read_tsv(PED,
                       col_names = "SAMPLE",
                       col_types = "i",
                       col_select = 1)

phenos = readr::read_csv(PHENO_FILE) %>% 
  # create sample column
  dplyr::mutate(SAMPLE = stringr::str_remove(fish, "KC"),
                SAMPLE = as.numeric(SAMPLE)) 

# Join

phen = pedp %>% 
  dplyr::left_join(.,
                   phenos,
                   by = "SAMPLE") %>% 
  dplyr::mutate(IID = SAMPLE) %>% 
  dplyr::select(SAMPLE, IID, all_of(PHENOTYPE))

readr::write_tsv(phen, OUT, col_names = F)
