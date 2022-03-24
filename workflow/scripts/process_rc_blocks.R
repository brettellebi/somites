# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
SITE_FILTER = "all_sites"
BIN_LENGTH = as.numeric("5000")
IN_FILE = file.path("/nfs/research/birney/users/ian/somites/recombination_blocks/F2",
                    SITE_FILTER,
                    paste(BIN_LENGTH, ".txt", sep = ""))

## True
IN_FILE = snakemake@input[[1]]
OUT_FILE = snakemake@output[[1]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()

# Read in HMM output file

df = readr::read_tsv(IN_FILE,
                     col_types = "ciiidii") %>% 
  # add key variables
  dplyr::mutate(SAMPLE = basename(sample) %>% 
                  stringr::str_remove(".txt") %>% 
                  as.numeric(.),
                CHROM = chr %>% 
                  as.numeric(),
                BIN_LENGTH = as.numeric(BIN_LENGTH),
                BIN_START = (bin - 1) * BIN_LENGTH + 1,
                BIN_END = bin * BIN_LENGTH,
                BIN_LENGTH_KB = BIN_LENGTH / 1e3,
                READS_PER_BIN = mat + pat) %>% 
  # Recode state to put in correct order (so that 0 == Cab)
  dplyr::mutate(STATE = dplyr::recode(state,
                                      `0` = 2,
                                      `1` = 1,
                                      `2` = 0)) %>% 
  # Rename columns
  dplyr::rename(CAB_READS = mat,
                KAGA_READS = pat)

# Read in total medaka genome count

## Get chromosome lengths
med_chr_lens = read.table("data/Oryzias_latipes.ASM223467v1.dna.toplevel.fa_chr_counts.txt",
                          col.names = c("chr", "end"))
## Add start
med_chr_lens$start = 1
## Reorder
med_chr_lens = med_chr_lens %>% 
  dplyr::select(chr, start, end) %>% 
  # remove MT
  dplyr::filter(chr != "MT") %>% 
  # convert chr to numeric
  dplyr::mutate(chr = as.numeric(chr))

# Fill in blocks based on genotype to the left

out = df %>% 
  # loop over SAMPLE
  split(., f = .$SAMPLE) %>% 
  purrr::map_dfr(., function(SAMPLE_DF){
    # loop over CHR
    SAMPLE_DF %>% 
      split(., f = .$CHROM) %>% 
      purrr::map_dfr(., function(CHR_DF){
        
        # Get target chromosome
        TARGET_CHROM = unique(CHR_DF$CHROM)
        # Get chromosome length
        CHR_LEN = med_chr_lens %>% 
          dplyr::filter(chr == TARGET_CHROM) %>% 
          dplyr::pull(end)
        # Create DF with all possible bins
        all_bins = tibble::tibble(CHROM = TARGET_CHROM,
                                  BIN = 1:ceiling( CHR_LEN / BIN_LENGTH),
                                  BIN_START = (BIN - 1) * BIN_LENGTH + 1,
                                  BIN_END = BIN * BIN_LENGTH)
        # Replace final bin end with CHR_LEN
        all_bins[nrow(all_bins), "BIN_END"] = CHR_LEN
        # Join ALL_BINS with df
        CHR_OUT = dplyr::left_join(all_bins,
                                   CHR_DF %>% 
                                     dplyr::select(SAMPLE,
                                                   CHROM,
                                                   BIN_START,
                                                   BIN_END,
                                                   CAB_READS,
                                                   KAGA_READS,
                                                   STATE),
                                   by = c("CHROM", "BIN_START", "BIN_END")) %>% 
          # fill SAMPLE column
          tidyr::fill(SAMPLE, .direction = "updown") %>% 
          # create new column with empty `STATE` values filled with previous value
          dplyr::mutate(STATE_IMP = STATE) %>% 
          tidyr::fill(STATE_IMP, .direction = "up")
        
        return(CHR_OUT)
        
      })
    
  })

# Write to file

readr::write_csv(out, OUT_FILE)
