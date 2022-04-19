# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#GENO_FILE = "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn_true/None/5000/1.csv"
#PHENOTYPE = "unsegmented_psm_area"

## True
GENO_FILE = snakemake@input[["genos"]]
OUT_PED = snakemake@output[["ped"]]
OUT_MAP= snakemake@output[["map"]]
OUT_PHEN = snakemake@output[["phen"]]
PHENOTYPE = snakemake@params[["phenotype"]]

# Read in genos

df = readr::read_csv(GENO_FILE)
  
# Impute missing genotypes based on previous call

df_imp = df %>% 
  dplyr::select(SAMPLE, CHROM, BIN, STATE) %>% 
  tidyr::pivot_wider(names_from = SAMPLE,
                     values_from = STATE) %>% 
  # group by chromosome so that values aren't pulled from previous chromosome
  dplyr::group_by(CHROM) %>% 
  # impute based on previous call. If first call is missing at start of chromosome, fill by genotypes below
  tidyr::fill(-c(CHROM, BIN), .direction = "updown") %>% 
  # ungroup
  dplyr::ungroup() %>% 
  # pivot back to longer format
  tidyr::pivot_longer(cols = -c(CHROM, BIN),
                      names_to = "SAMPLE",
                      values_to = "STATE") %>% 
  # make `SAMPLE` numeric
  dplyr::mutate(SAMPLE = as.numeric(SAMPLE))

# Process

new = df_imp %>%
  #dplyr::slice_sample(n = 1000) %>% 
  # order sample
  dplyr::arrange(SAMPLE, CHROM, BIN) %>% 
  # add SNP ID by pasting CHROM and BIN_START
  dplyr::mutate(SNP = paste(CHROM, BIN, sep = ":")) %>% 
  dplyr::select(SNP, CHROM, BIN, SAMPLE, STATE) %>% 
  # split genotype into two columns, one for each allele
  dplyr::mutate(GT = dplyr::case_when(STATE == 0 ~ "AA",
                                      STATE == 1 ~ "AB",
                                      STATE == 2 ~ "BB"))

# Create .ped file by widening to put genotypes for each SNP into their own columns

ped = new %>% 
  dplyr::select(SAMPLE, SNP, GT) %>% 
  # put each sample's genotypes into their own column
  tidyr::pivot_wider(names_from = SNP, values_from = GT, names_vary = "slowest")

# Create .map file

map = new %>% 
  # order
  dplyr::arrange(CHROM, BIN) %>% 
  # get distinct loci
  dplyr::distinct(CHROM, SNP, BIN, .keep_all = F) %>% 
  # reorder columns
  dplyr::select(CHROM, SNP, BIN)

# Create .phen file (see example here: https://github.com/jianyangqt/gcta/blob/master/test/tests/data/test.phen )

#phen = new %>% 
#  dplyr::distinct(FID, IID, all_of(PHENOTYPE)) %>% 
#  dplyr::arrange(FID)

# Write to file

readr::write_tsv(ped, OUT_PED, col_names = F)
readr::write_tsv(map, OUT_MAP, col_names = F)
#readr::write_tsv(phen, OUT_PHEN, col_names = F)
