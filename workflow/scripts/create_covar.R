# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#PED = "/hps/nobackup/birney/users/ian/somites/peds/F2/hdrr/None/5000/1/intercept.ped"
#PHENO_FILE = here::here("data/20220321_phenotypes.xlsx")
#COVARS = "Microscope"

## True
PED = snakemake@input[["ped"]]
PHENO_FILE = snakemake@input[["phenos"]]
COVARS = snakemake@params[["covars"]]
OUT = snakemake@output[[1]]

# Read in files

pedp = readr::read_tsv(PED,
                       col_names = "SAMPLE",
                       col_types = "i",
                       col_select = 1)

phenos = readr::read_delim(PHENO_FILE, delim =";")

# Process covariates string

covs = stringr::str_split(COVARS, "-") %>% 
  unlist()

# Join

covar = dplyr::left_join(pedp,
                         phenos,
                         by = "SAMPLE") %>% 
  dplyr::mutate(IID = SAMPLE) %>% 
  dplyr::select(SAMPLE, IID, all_of(covs))

readr::write_tsv(covar, OUT, col_names = F)
