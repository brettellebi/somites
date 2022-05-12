# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug

GENOS = "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn_true/None/5000/0.8.csv"
PHENOS = "/hps/nobackup/birney/users/ian/somites/phens/inv_norm/intercept.phen"
RES = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/true/hdrr/None/5000/0.8/intercept.loco.mlma"
MIN_P = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/min_p/hdrr/None/5000/0.8/intercept.csv"
BIN_LENGTH = 5000

## True
GENOS = snakemake@input[["genos"]]
PHENOS = snakemake@input[["phenos"]]
RES = snakemake@input[["res"]]
MIN_P = snakemake@input[["min_p"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT = snakemake@output[["fig"]]

######################
# Read in files
######################

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

# Genos

genos_raw = readr::read_csv(GENOS) %>% 
  # add bin start and end
  dplyr::mutate(BIN_START = (BIN * BIN_LENGTH) + 1,
                BIN_END = (BIN + 1) * BIN_LENGTH) 

genos = genos_raw %>% 
  # filter for significant region
  dplyr::filter(CHROM == 3 & BIN_START >= loc_start & BIN_END <= loc_end) %>% 
  dplyr::mutate(BIN_START_MB = BIN_START / 1e6 )

# Phenos

phenos = readr::read_tsv(PHENOS,
                         col_names = c("FID", "SAMPLE", "split_invnorm_intercept")) %>% 
  dplyr::select(-FID)


######################
# Get mean proportion of Kaga per sample across locus
######################

df = genos %>% 
  dplyr::group_by(SAMPLE) %>% 
  dplyr::summarise(MEAN_PROP_KAGA = mean(PROP_KAGA)) %>% 
  # bind with phenos
  dplyr::left_join(phenos,
                   by = "SAMPLE")

######################
# Plot
######################

genos_raw %>% 
  dplyr::slice_sample(n = 1e6) %>% 
  dplyr::mutate(STATE = factor(STATE, levels = c(0,1,2))) %>% 
  ggplot() +
  geom_boxplot(aes(STATE, PROP_KAGA, group = STATE)) +
  theme_bw()

out = df %>% 
  ggplot(aes(MEAN_PROP_KAGA, split_invnorm_intercept)) +
  geom_point() +
  geom_smooth(method = "lm") +
  cowplot::theme_cowplot() +
  ggtitle("Chr3 significant region:\nMean proportion of Kaga reads across locus\nvs split inverse-normalised period intercept") +
  xlab("mean proportion of Kaga reads across locus") +
  ylab("split inverse-normalised period intercept")


ggsave(OUT,
       out,
       device = "png",
       width = 10,
       height = 5,
       units = "in", 
       dpi = 400)
