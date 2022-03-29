# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#GENO_FILE = "/hps/nobackup/birney/users/ian/somites/processed_recomb/F2/all_sites/20000.csv"
#PHENO_FILE = here::here("data/20220321_phenotypes.xlsx")
#PHENOTYPE = "unsegmented_psm_area"

## True
GENO_FILE = snakemake@input[["genos"]]
PHENO_FILE = snakemake@input[["phenos"]]
OUT_PED = snakemake@output[["ped"]]
OUT_MAP= snakemake@output[["map"]]
OUT_PHEN = snakemake@output[["phen"]]
PHENOTYPE = snakemake@params[["phenotype"]]

# Read in files

df = readr::read_csv(GENO_FILE)

phenos = readxl::read_xlsx(PHENO_FILE) %>% 
  dplyr::select(SAMPLE, all_of(PHENOTYPE))

# Process

new = df %>%
  #dplyr::slice_sample(n = 1000) %>% 
  # order sample
  dplyr::arrange(SAMPLE) %>% 
  dplyr::mutate(SNP = paste(CHROM, BIN_START, sep = ":")) %>% 
  dplyr::select(SNP, CHROM, BIN_START, SAMPLE, STATE_IMP) %>% 
  # add extra columns
  dplyr::mutate(FID = SAMPLE,
                IID = SAMPLE,
                IID_PAT = 0,   # dummy
                IID_MAT = 0,   # dummy
                SEX = 0        # dummy
                ) %>% 
  # bind with phenotypes
  dplyr::left_join(.,
                   phenos,
                   by = "SAMPLE") %>% 

  # split genotype into two columns, one for each allele
  dplyr::mutate(A1 = dplyr::case_when(STATE_IMP == 0 ~ "A",
                                      STATE_IMP == 1 ~ "A",
                                      STATE_IMP == 2 ~ "B"),
                A2 = dplyr::case_when(STATE_IMP == 0 ~ "A",
                                      STATE_IMP == 1 ~ "B",
                                      STATE_IMP == 2 ~ "B"))

# Create .ped file by widening to put genotypes for each SNP into their own columns

ped = new %>% 
  dplyr::select(FID, IID, IID_PAT, IID_MAT, SEX, all_of(PHENOTYPE), SNP, A1, A2) %>% 
  # put each sample's genotypes into their own column
  tidyr::pivot_wider(names_from = SNP, values_from = c(A1, A2), names_vary = "slowest")

# Create .map file

map = new %>% 
  dplyr::select(CHROM, SNP, BIN_START) %>% 
  # add position in centimorgans column with dummy value of 0
  dplyr::mutate(POS = 0) %>% 
  # reorder
  dplyr::select(CHROM, SNP, POS, BIN_START)

# Create .phen file (see example here: https://github.com/jianyangqt/gcta/blob/master/test/tests/data/test.phen)

phen = new %>% 
  dplyr::select(FID, IID, all_of(PHENOTYPE))

# Write to file

readr::write_delim(ped, OUT_PED, delim = " ", col_names = F)
readr::write_delim(map, OUT_MAP, delim = " ", col_names = F)
readr::write_delim(phen, OUT_PHEN, delim = " ", col_names = F)
  
