# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(circlize)

# Get variables

## Debug
#IN_FILE = "/hps/nobackup/birney/users/ian/somites/genos/F0_and_F1/counts/Cab/5000.csv"
#CHROM_LENGTHS = here::here("config/hdrr_chrom_lengths.csv")
#SAMPLE = "Cab"
#BIN_LENGTH = 5000
#PAL = "#43AA8B"

## True
IN_FILE = snakemake@input[["gt_counts"]]
CHROM_LENGTHS = snakemake@input[["chrom_lens"]]
SAMPLE = snakemake@params[["sample"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
PAL = snakemake@params[["palette"]]
OUT_PLOT = snakemake@output[["plot"]]

# Get lighter/darker functions

source("https://gist.githubusercontent.com/brettellebi/c5015ee666cdf8d9f7e25fa3c8063c99/raw/91e601f82da6c614b4983d8afc4ef399fa58ed4b/karyoploteR_lighter_darker.R")


# Read in data

genos = readr::read_csv(IN_FILE)

chroms = readr::read_csv(CHROM_LENGTHS,
                         col_names = c("CHROM", "LENGTH")) %>% 
  # remove MT
  dplyr::filter(CHROM != "MT") %>% 
  dplyr::mutate(CHROM = CHROM %>% 
                  stringr::str_replace(pattern = "^", replacement = "chr"),
                START = 1) %>% 
  dplyr::select(CHROM, START, END = LENGTH)


# Process `genos`

if (SAMPLE == "F1"){
  TARGET_CHROM = "PROP_HET"
} else {
  TARGET_CHROM = "PROP_HOM"
}

homozyg = genos %>% 
  dplyr::mutate(CHROM = CHROM %>% 
                  stringr::str_replace(pattern = "^", replacement = "chr"),
                BIN_START = (BIN -1) * BIN_LENGTH + 1,
                BIN_END = BIN * BIN_LENGTH) %>% 
  dplyr::select(CHROM, BIN_START, BIN_END, dplyr::all_of(TARGET_CHROM))

n_vars = genos %>% 
  dplyr::mutate(CHROM = CHROM %>% 
                  stringr::str_replace(pattern = "^", replacement = "chr"),
                BIN_START = (BIN -1) * BIN_LENGTH + 1,
                BIN_END = BIN * BIN_LENGTH) %>% 
  dplyr::select(CHROM, BIN_START, BIN_END, TOT_HITS) 

# Set output

png(OUT_PLOT,
    width = 20,
    height = 20,
    units = "cm",
    res = 500)

# Create Circos plots

circos.par(cell.padding = c(0, 0, 0, 0),
           track.margin = c(0, 0),
           gap.degree = c(rep(1, nrow(chroms) - 1), 8))

# Initialize plot

circos.initializeWithIdeogram(chroms,
                              plotType = c("axis", "labels"),
                              major.by = 1e7,
                              axis.labels.cex = 0.25*par("cex"))

if (SAMPLE == "F1"){
  CENTER_LAB = paste(SAMPLE,
                     "\nheterozygosity\nand\nvariant count",
                     "\nwithin\n",
                     BIN_LENGTH / 1000,
                     "kb bins",
                     sep = "")
} else {
  CENTER_LAB = paste(SAMPLE,
                     "\nhomozygosity\nand\nvariant count",
                     "\nwithin\n",
                     BIN_LENGTH / 1000,
                     "kb bins",
                     sep = "")
}
# Add label to center
text(0, 0, CENTER_LAB)

# Add proportion of homozygosity

circos.genomicTrack(homozyg,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          type = "h",
                                          col = PAL,
                                          cex = 0.05)
                    },
                    track.height = 0.12,
                    bg.border = NA,
                    ylim = c(0, 1))
# y-axis label
circos.yaxis(side = "right",
             at = c(0, 0.5, 1),
             labels.cex = 0.25*par("cex"),
             tick.length = 2
)
# y-axis title

if (SAMPLE == "F1"){
  AXIS_LAB = "proportion\nheterozygous"
} else {
  AXIS_LAB = "proportion\nhomozygous"
}

circos.text(0, 0.25,
            labels = AXIS_LAB,
            sector.index = "chr1",
            facing = "clockwise",
            adj = c(0, -0.5),
            cex = 0.3*par("cex"))

# Add number of hits

## get max number of variants

MAX_VARS = max(n_vars$TOT_HITS, na.rm = T)

circos.genomicTrack(n_vars,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          type = "h",
                                          col = "#F3B700",
                                          cex = 0.05)
                    },
                    track.height = 0.12,
                    bg.border = NA,
                    ylim = c(0, MAX_VARS))
# y-axis label
circos.yaxis(side = "right",
             at = c(0, 500),
             labels.cex = 0.25*par("cex"),
             tick.length = 2
)

circos.text(0, 0,
            labels = "N variants\nper bin",
            sector.index = "chr1",
            facing = "clockwise",
            adj = c(0, -0.5),
            cex = 0.3*par("cex"))
