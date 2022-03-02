# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
GWAS_RESULTS = "/hps/nobackup/birney/users/ian/somites/association_testing/20220214/all_sites/true_results/unsegmented_psm_area/None/FALSE/5000.rds"
SIG_LEVELS = "data/20220214_permutation_mins.csv" # True phenotypes
SOURCE_FILE = "workflow/scripts/get_manhattan_source.R"
SITE_FILTER = "all_sites"
BIN_LENGTH = 5000
TARGET_PHENO = "unsegmented_psm_area"
COVARIATES = "None"
INVERSE_NORM = "FALSE"

## True
### Input
GWAS_RESULTS = snakemake@input[["gwas_results"]]
SIG_LEVELS = snakemake@input[["sig_levels"]]
SOURCE_FILE = snakemake@input[["source_file"]]
### Parameters
SITE_FILTER = snakemake@params[["site_filter"]]
TARGET_PHENO = snakemake@params[["target_phenotype"]]
COVARIATES = snakemake@params[["covariates"]]
INVERSE_NORM = snakemake@params["inverse_norm"]
BIN_LENGTH = snakemake@params[["bin_length"]] %>%
  as.numeric()
### Output
OUT_FILE = snakemake@output[["fig"]]

# Read in source file

source(SOURCE_FILE)

# Read in GWAS results

RESULTS = readRDS(GWAS_RESULTS)

# Read in permutation minimums

SIG_LEVEL = readr::read_csv(SIG_LEVELS,
                            col_types = c("ccccid")) %>% 
  dplyr::filter(SITE_FILTER = SITE_FILTER,
                TARGET_PHENO = TARGET_PHENO,
                COVARIATES = COVARIATES,
                INVERSE_NORM = INVERSE_NORM,
                BIN_LENGTH = BIN_LENGTH) %>% 
  dplyr::pull(MIN_P)

# Choose palette

pal = eval(as.name(paste(TARGET_PHENO, "_pal", sep = "")))

# Generate Manhattan plot

out_clean = clean_gwas_res(RESULTS,
                           bin_length = BIN_LENGTH,
                           chr_lens = med_chr_lens)

# Plot
out_plot = plot_man(out_clean,
                    site_filter = "all_sites",
                    phenotype = TARGET_PHENO,
                    bin_length = BIN_LENGTH, 
                    gwas_pal = pal,
                    med_chr_lens = med_chr_lens,
                    sig_level = SIG_LEVEL) +
  ylim(0,7) + 
  labs(subtitle = paste("Covariates: ", COVARIATES,
                        "\nInverse normalised: ", INVERSE_NORM,
                        "\nn samples: ", N_SAMPLES,
                        sep = ""))

# Save

ggsave(OUT_FILE,
       out_plot,
       device = "png",
       width = 9.6,
       height = 6,
       units = "in",
       dpi = 400)
