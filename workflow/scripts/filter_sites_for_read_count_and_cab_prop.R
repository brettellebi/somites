#IN_DIR = "/hps/nobackup/birney/users/ian/somites/dpABs/F2/no_repeat_reads_or_pers_hets"
#SITES_FILE = "/hps/nobackup/birney/users/ian/somites/data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats_no_persistent_hets.txt"
#IN_FILES = list.files(IN_DIR, full.names = T)[1:4]
#MIN_READS = 100
#MAX_READS = 10000
#MIN_PROP_CAB = 0.2
#MAX_PROP_CAB = 0.8

IN_FILES = snakemake@input[["dp_files"]]
SITES_FILE = snakemake@input[["sites_file"]]
OUT_FILE_COUNTS_PRE = snakemake@output[["counts_pre_filter"]]
OUT_FILE_COUNTS_POST = snakemake@output[["counts_post_filter"]]
OUT_FILE_SITES = snakemake@output[["filtered_sites"]]
LOG_FILE = snakemake@log[[1]]
MIN_READS = snakemake@params[["min_reads"]]
MAX_READS = snakemake@params[["max_reads"]]
MIN_PROP_CAB = snakemake@params[["min_prop_cab"]]
MAX_PROP_CAB = snakemake@params[["max_prop_cab"]]

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

# Group by site and get total CAB and KAGA reads
counts_df = df %>% 
  dplyr::select(-c(CAB_ALLELE, KAGA_ALLELE)) %>% 
  dplyr::group_by(CHROM, POS) %>% 
  dplyr::summarise(TOTAL_CAB = sum(CAB),
                   TOTAL_KAGA = sum(KAGA)) %>% 
  dplyr::ungroup() %>% 
  # get overall total count and proportion of CAB
  dplyr::mutate(TOTAL_READS = TOTAL_CAB + TOTAL_KAGA,
                PROP_CAB = TOTAL_CAB / TOTAL_READS) 

# Write pre-filter counts and proportions to file
readr::write_tsv(counts_df, OUT_FILE_COUNTS_PRE)

counts_df = counts_df %>% 
  # apply filters
  ## total reads
  dplyr::filter(TOTAL_READS >= MIN_READS & TOTAL_READS <= MAX_READS) %>% 
  ## proportion Cab
  dplyr::filter(PROP_CAB >= MIN_PROP_CAB & PROP_CAB <= MAX_PROP_CAB) %>% 
  # add CHR:POS column
  tidyr::unite(col = "CHROM:POS", c(CHROM, POS), sep = ":", remove = F)

# Write post-filter counts and proportions to file
counts_df %>% 
  dplyr::select(-`CHROM:POS`) %>% 
  readr::write_tsv(OUT_FILE_COUNTS_POST)

# Read in sites file, filter, and write to file
sites_df = readr::read_tsv(SITES_FILE,
                           col_names = c("CHROM", "POS", "POS_2", "REF", "ALT", "CAB_GT", "KAGA_GT"),
                           col_types = c("iiicccc")) %>% 
  
  # add CHR:POS column
  tidyr::unite(col = "CHROM:POS", c(CHROM, POS), sep = ":", remove = F) %>% 
  # filter for sites in `counts_df`
  dplyr::filter(`CHROM:POS` %in% counts_df$`CHROM:POS`) %>% 
  # remove `CHROM:POS` column
  dplyr::select(-`CHROM:POS`) %>% 
  # write to file
  readr::write_tsv(OUT_FILE_SITES, col_names = F)

