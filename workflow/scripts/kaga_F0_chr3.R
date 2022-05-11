# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug

GENOS = "/hps/nobackup/birney/users/ian/somites/genos/F0_and_F1/hdrr/final/all.csv"
COUNTS = "/hps/nobackup/birney/users/ian/somites/genos/F0_and_F1/hdrr/counts/Kaga/5000.csv"
RES = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/true/hdrr/None/5000/0.8/intercept.loco.mlma"
MIN_P = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/min_p/hdrr/None/5000/0.8/intercept.csv"

BIN_LENGTH = 5000

## True


######################
# Read in files
######################

# Genotypes

#genos = readr::read_csv(GENOS)

# Select columns

#genos = genos %>% 
#  dplyr::select(CHROM = `# [1]CHROM`,
#                POS = `[2]POS`,
#                KAGA_GT = `[7]Kaga:GT`)
#
# Counts

counts = readr::read_csv(COUNTS) %>% 
  dplyr::mutate(BIN_START = (BIN -1) * BIN_LENGTH + 1,
                BIN_END = BIN * BIN_LENGTH) %>% 
  # re-calculate PROP_HOM
  dplyr::mutate(dplyr::across(c("HOM", "HET", "MISS"),
                ~ tidyr::replace_na(.x, 0))) %>% 
  dplyr::mutate(TOT_HITS = HOM + HET + MISS) %>% 
  dplyr::mutate(PROP_HOM = HOM/ TOT_HITS)

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

# Filter genotypes

#genos_sig = genos %>% 
#  dplyr::filter(CHROM == 3 & POS >= loc_start & POS <= loc_end) %>% 
#  # Remove missing genotypes
#  dplyr::filter(!(KAGA_GT %in% c("./.", ".|."))) %>% 
#  # Recode GTs as hom or het
#  dplyr::mutate(HOM_HET = dplyr::case_when(KAGA_GT %in% c("0/1", "0|1") ~ "HET",
#                                           KAGA_GT %in% c("1|1", "0/0", "0|0", "1/1") ~ "HOM")) %>% 
#  # convert POS to Mb
#  dplyr::mutate(POS_MB = POS/1e6)

# Filter counts

counts_sig = counts %>% 
  dplyr::filter(CHROM == 3 & BIN_START >= loc_start & BIN_END <= loc_end) %>% 
  dplyr::mutate(BIN_START_MB = BIN_START / 1e6 )

######################
# Plot
######################

out = counts_sig %>% 
  ggplot() + 
  #geom_line(aes(BIN_START_MB, PROP_HOM)) +
  geom_area(aes(BIN_START_MB, PROP_HOM), fill = "#DE3C4B") +
  cowplot::theme_cowplot() +
  ggtitle("Chr3 significant region:\nKaga F0 proportion of homozygous SNPs within 5 kb bins") +
  xlab("position (Mb)") +
  ylab("proportion of homozygous SNPs")

ggsave(OUT,
       out,
       device = "png",
       width = 9,
       height = 3.2,
       units = "in", 
       dpi = 400)
