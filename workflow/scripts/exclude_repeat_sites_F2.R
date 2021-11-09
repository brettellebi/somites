# Fetch variables

#HOMO_DIV_FILE = "/hps/nobackup/birney/users/ian/somites/data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"
#REPEATS_FILE = "/nfs/research/birney/users/ian/mikk_genome/repeats/medaka_hdrr_repeats.fixed.gff"
TARGETS_FILE = snakemake@input[["target_sites"]]
REPEATS_FILE = snakemake@input[["repeats_file"]]
OUT_FILE = snakemake@output[[1]]
LOG_FILE = snakemake@log[[1]]

# Send log

log <- file(LOG_FILE, open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(GenomicRanges)

# Read in data

## Homozygous-divergent sites

target_df = readr::read_delim(TARGETS_FILE,
                              delim = "\t",
                              col_names = c("CHROM", "POS_1", "POS_2", "REF", "ALT", "GT_1", "GT_2"),
                              col_types = c("iiicccc"))

## Repeats file

hdrr_reps = read.table(REPEATS_FILE,
                       header = F, sep = "\t", skip = 3, comment.char = "", quote = "", as.is = T) %>%
  dplyr::select(chr = V1, start = V4, end = V5)

# Convert to GRanges

target_ranges = GenomicRanges::makeGRangesFromDataFrame(target_df,
                                                        keep.extra.columns = T,
                                                        seqnames.field = "CHROM",
                                                        start.field = "POS_1",
                                                        end.field = "POS_2")

repeat_ranges = GenomicRanges::makeGRangesFromDataFrame(hdrr_reps,
                                                        keep.extra.columns = T,
                                                        seqnames.field = "chr",
                                                        start.field = "start",
                                                        end.field = "end")

# Find overlaps

olaps_out = GenomicRanges::findOverlaps(target_ranges, repeat_ranges)

# How many sites overlap repeat regions? 
print(paste("NUMBER OF SITES OVERLAPPING REPEAT REGIONS:",
            length(unique(olaps_out@from))))

# What proportion of target regions are to be excluded?
print(paste("PROPORTION OF TARGET REGIONS TO BE EXCLUDED:",
            length(unique(olaps_out@from)) / nrow(target_df)))

# Filter for sites that don't overlap repeats 

target_out = target_df[-olaps_out@from, ]

# Write to file

readr::write_delim(target_out, OUT_FILE, delim = "\t", col_names = F)
