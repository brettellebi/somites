# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#GWAS_RESULTS = "/hps/nobackup/birney/users/ian/somites/association_testing/20220214/all_sites/true_results/intercept/Microscope/TRUE/20000.rds"
#PERM_RESULTS = list("/hps/nobackup/birney/users/ian/somites/association_testing/20220214/all_sites/permutations/intercept/Microscope/TRUE/20000/2.rds",
#                    "/hps/nobackup/birney/users/ian/somites/association_testing/20220214/all_sites/permutations/intercept/Microscope/TRUE/20000/8.rds")
#SOURCE_FILE = "workflow/scripts/get_manhattan_source.R"
#SITE_FILTER = "all_sites"
#BIN_LENGTH = 20000
#TARGET_PHENO = "intercept"
#COVARIATES = "Microscope"
#INVERSE_NORM = "TRUE"

## True
### Input
GWAS_RESULTS = snakemake@input[["gwas_results"]]
PERM_RESULTS = snakemake@input[["perm_results"]]
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

# Get signficance level for bonferroni correction

ALPHA = 0.05

# Read in source file

source(SOURCE_FILE)

# Read in GWAS results

RESULTS = readRDS(GWAS_RESULTS)

# Get relevant p-value column (if no covariates are specified, it's `p_value_REML`)
# Otherwise it's `p_value_REML.1`

if ("p_value_REML" %in% colnames(RESULTS$results)){
  P_COL = "p_value_REML"
} else {
  P_COL = "p_value_REML.1"
}

# Rename column in results

if (P_COL == "p_value_REML.1"){
  RESULTS$results = RESULTS$results %>% 
    dplyr::rename(p_value_REML = p_value_REML.1)
}

# Read in permutation results and get `SIG_LEVEL`

PERM_LIST = purrr::map(PERM_RESULTS, function(PERM){
  readRDS(PERM)
})
names(PERM_LIST) = PERM_RESULTS %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".rds")

perm_df = purrr::map_dfr(PERM_LIST, function(PERM){
  OUT = tibble::tibble(MIN_P = PERM$results %>% 
                         dplyr::select(dplyr::all_of(P_COL)) %>%
                         min(., na.rm = T)
  )
}, .id = "SEED")

# Get minimum
SIG_LEVEL = min(perm_df$MIN_P)

# Choose palette

pal = eval(as.name(paste(TARGET_PHENO, "_pal", sep = "")))

# Generate Manhattan plot

## Prepare data
out_clean = clean_gwas_res(RESULTS,
                           bin_length = BIN_LENGTH,
                           chr_lens = med_chr_lens)

## Get bonferroni significance level
BONFERRONI = ALPHA / nrow(out_clean)

# Plot
out_plot = plot_man(out_clean,
                    site_filter = SITE_FILTER,
                    phenotype = TARGET_PHENO,
                    bin_length = BIN_LENGTH, 
                    gwas_pal = pal,
                    med_chr_lens = med_chr_lens,
                    sig_level = SIG_LEVEL,
                    bonferroni = BONFERRONI) +
  labs(subtitle = paste("Covariates: ", COVARIATES,
                        "\nInverse-normalised: ", INVERSE_NORM,
                        sep = ""))

# Save
## Make sure the directory exists
dir.create(dirname(OUT_FILE), recursive = T, showWarnings = F)
## Write
ggsave(OUT_FILE,
       out_plot,
       device = "png",
       width = 9.6,
       height = 6,
       units = "in",
       dpi = 400)
