# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#GENO_FILE = "/hps/nobackup/birney/users/ian/somites/genos/F0_and_F1/hdrr/final/all.csv"
#CHROM_LENGTHS = here::here("config/hdrr_chrom_lengths.csv")
#SAMPLE = "Cab"
#BIN_LENGTH = 5000

## True
GENO_FILE = snakemake@input[["genos"]]
CHROM_LENGTHS = snakemake@input[["chrom_lens"]]
SAMPLE = snakemake@params[["sample"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT_FILE = snakemake@output[[1]]

# Read in files

chr_lens = readr::read_csv(CHROM_LENGTHS,
                           col_names = c("CHROM", "END")) %>% 
  # remove MT
  dplyr::filter(CHROM != "MT") %>% 
  dplyr::mutate(CHROM = as.numeric(CHROM))

genos = readr::read_csv(GENO_FILE)

# Tidy up column names
colnames(genos) = colnames(genos) %>% 
  stringr::str_remove(pattern = "\\[.\\]") %>% 
  stringr::str_remove(pattern = ":GT") %>% 
  stringr::str_remove(pattern = "# ")

# Process

out = genos %>% 
  # rename sample column
  dplyr::select(CHROM, POS, GT = all_of(SAMPLE)) %>% 
  # remove dash and line
  dplyr::mutate(GT_NEW = dplyr::case_when(GT == "1|1" | GT == "1/1" ~ "11",
                                          GT == "0|0" | GT == "0/0" ~ "00",
                                          GT == "0|1" | GT == "0/1" ~ "01",
                                          GT == ".|." | GT == "./." ~ "MISS",
                                          TRUE ~ GT),
                HOM_HET = dplyr::case_when(GT_NEW == "11" | GT_NEW == "00" ~ "HOM",
                                           GT_NEW == "01" ~ "HET",
                                           GT_NEW == "MISS" ~ GT_NEW,
                                           TRUE ~ GT_NEW))


# Get bin labels
df_bins = purrr::map_dfr(chr_lens$CHROM, function(TARGET_CHROM){
  # Filter DF for that chromosome
  df_filt = out %>% 
    dplyr::filter(CHROM == TARGET_CHROM)
  
  # Get length of target chromosome
  TARGET_END = chr_lens %>% 
    dplyr::filter(CHROM == TARGET_CHROM) %>% 
    dplyr::pull(END)
  
  # Create sequence of breaks
  BINS = seq(1, TARGET_END, by = BIN_LENGTH)
  BINS = c(BINS, TARGET_END)
  
  # Create column with bins
  result = df_filt %>% 
    dplyr::mutate(BIN = cut(POS, breaks = BINS, right = F, labels = F))
  
  return(result)
})


# Group by chromosome and bin, and get total variants within each bin

final = df_bins %>% 
  dplyr::count(CHROM, BIN, HOM_HET) %>% 
  tidyr::pivot_wider(names_from = "HOM_HET",
                     values_from = n) %>% 
  # replace NAs with 0
  dplyr::mutate(dplyr::across(c("MISS", "HOM", "HET"),
                              ~tidyr::replace_na(.x, 0))) %>% 
  dplyr::mutate(TOT_HITS = HOM + HET,
                PROP_HOM = HOM / TOT_HITS,
                PROP_HET = HET / TOT_HITS) %>% 
  # remove rows with NA (caused by bins with only missing genotypes)
  tidyr::drop_na()

# Save to file

readr::write_csv(final, OUT_FILE)
