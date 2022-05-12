# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug

COUNTS = "/hps/nobackup/birney/users/ian/somites/genos/F0_and_F1/hdrr/counts/Kaga/5000.csv"
RES = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/true/hdrr/None/5000/0.8/intercept.loco.mlma"
MIN_P = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/min_p/hdrr/None/5000/0.8/intercept.csv"
BIN_LENGTH = 5000

## True
COUNTS = snakemake@input[["counts"]]
RES = snakemake@input[["res"]]
MIN_P = snakemake@input[["min_p"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT = snakemake@output[["fig"]]

######################
# Read in files
######################

# Counts

counts = readr::read_csv(COUNTS) %>% 
  dplyr::mutate(BIN_START = (BIN -1) * BIN_LENGTH + 1,
                BIN_END = BIN * BIN_LENGTH)

# Significance level

PERM_SIG = readr::read_csv(MIN_P) %>% 
  dplyr::pull(MIN_P) %>% 
  min(.)

# GWAS results

res = readr::read_tsv(RES) %>% 
  # Filter for significant bins
  dplyr::filter(p < PERM_SIG) %>% 
  # Filter for chr3
  dplyr::filter(Chr == 3) %>% 
  # Add POS
  dplyr::mutate(BIN_START = (bp * BIN_LENGTH) + 1,
                BIN_END = (bp + 1) * BIN_LENGTH)

# Get range

loc_start = min(res$BIN_START)
loc_end = max(res$BIN_END)

# Filter counts

counts_sig = counts %>% 
  dplyr::filter(CHROM == 3 & BIN_START >= loc_start & BIN_END <= loc_end) %>% 
  dplyr::mutate(BIN_START_MB = BIN_START / 1e6 )

######################
# Plot
######################

prop_hom = counts_sig %>% 
  ggplot() + 
  #geom_line(aes(BIN_START_MB, PROP_HOM)) +
  geom_area(aes(BIN_START_MB, PROP_HOM), fill = "#DE3C4B") +
  cowplot::theme_cowplot() +
  ggtitle("Chr3 significant region:\nKaga F0 proportion of homozygous SNPs within 5 kb bins") +
  xlab("position (Mb)") +
  ylab("proportion of homozygous SNPs")

snp_counts = counts_sig %>% 
  ggplot() +
  geom_col(aes(BIN_START_MB, TOT_HITS), fill = "#F3B700") + 
  cowplot::theme_cowplot() +
  ggtitle("SNP counts per bin") +
  xlab("position (Mb)") +
  ylab("SNP count per bin")

out = cowplot::plot_grid(prop_hom,
                         snp_counts,
                         rel_widths = 1,
                         rel_heights = 0.5,
                         nrow = 2)

ggsave(OUT,
       out,
       device = "png",
       width = 10,
       height = 5,
       units = "in", 
       dpi = 400)
