# Fetch variables

#IN_SITES =  "/hps/nobackup/birney/users/ian/somites/data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats.txt"
#BLACK_LIST = "/hps/software/users/birney/ian/repos/somites/data/20200716_mikk_panel_proportion_het_across_genome.txt"

IN_SITES = snakemake@input[["sites_excl_repeats"]]
BLACK_LIST = snakemake@input[["black_list"]]
OUT_FILE = snakemake@output[[1]]
LOG_FILE = snakemake@log[[1]]

# Send log

log <- file(LOG_FILE, open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(GenomicRanges)

# Read in homozygous-divergent sites file

homdiv = readr::read_tsv(IN_SITES,
                         col_names = c("CHROM", "START", "END", "F0_1_ALLELE", "F0_2_ALLELE", "F0_GT", "F1_GT"),
                         col_types = c("iiicccc"))

# Read in black list

black = readr::read_tsv(BLACK_LIST,
                        col_types = c("iiid")) %>% 
  # filter for blocks with > 90% heterozygosity
  dplyr::filter(p_het > 0.9)

# Convert both to ranges

## hom-div sites
homdiv_ranges = homdiv %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM",
                                          start.field = "START",
                                          end.field = "END")

## black sites
black_ranges = black %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "chr",
                                          start.field = "start",
                                          end.field = "end")

# Find overlaps

olaps = GenomicRanges::findOverlaps(homdiv_ranges,
                                    black_ranges)

# Filter out hom-div sites which overlap black list

out = homdiv %>% 
  dplyr::slice(-olaps@from)

# Save file

readr::write_tsv(out, OUT_FILE, col_names = F)
