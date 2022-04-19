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

## True

PHEN = snakemake@input[["phen"]]
COVAR = snakemake@input[["covar"]]
OUT = snakemake@output[[1]]

# Read in files and combine

phen = readr::read_tsv(PHEN,
                      col_names = c("FID", "IID", "PHEN"))

covar = readr::read_tsv(COVAR,
                        col_names = c("FID", "IID", "Microscope"))

df = left_join(phen, covar, by = c("FID", "IID"))
  
# Create function
  
my_invnorm = function(x) {
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}
  
# Split by microscope and inverse-normalise

phen_au = df %>% 
  dplyr::filter(Microscope == "AU") %>% 
  dplyr::mutate(INV_NORM = my_invnorm(PHEN))

phen_db = df %>% 
  dplyr::filter(Microscope == "DB") %>% 
  dplyr::mutate(INV_NORM = my_invnorm(PHEN))

# Combine

invnorm = dplyr::bind_rows(phen_au,
                           phen_db)

# Bind with original phen to ensure the same order

out = phen %>% 
  dplyr::left_join(.,
                   invnorm,
                   by = c("FID", "IID")) %>% 
  # select and rename columns
  dplyr::select(FID, IID, PHEN = INV_NORM)

# Write to file

readr::write_tsv(out, OUT, col_names = F)
