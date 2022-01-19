# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(GridLMM)
library(KRLS)

# Get variables

## Debugging
GENO_FILE = "/nfs/research/birney/users/ian/somites/association_testing/20220118/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/inputs/20000.rds"
PHENO_FILE = "data/20210917_First400_F2_DF.xlsx"
SOURCE_FILE = "workflow/scripts/run_gwls_source.R"
BIN_LENGTH = 20000
TARGET_PHENO = "intercept"
PERM_SEED = 1

## True
GENO_FILE = snakemake@input[["gt_pos_list"]]
PHENO_FILE = snakemake@input[["phenotypes_file"]]
SOURCE_FILE = snakemake@input[["source_file"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>%
    as.numeric()
TARGET_PHENO = snakemake@params[["target_phenotype"]]
PERM_SEED = snakemake@params[["permutation_seed"]] %>%
    as.numeric()
OUT_FILE = snakemake@output[[1]]

# Get GWAS functions

source(SOURCE_FILE)

# Load genotypes and positions

in_list = readRDS(GENO_FILE)

# Read in phenotypes

## Set seed for randomisation of phenotype

set.seed(PERM_SEED)

## Read in file and wrangle
phenos = readxl::read_xlsx(PHENO_FILE) %>%
    # adjust sample names
    dplyr::mutate(SAMPLE = fish %>% stringr::str_remove("KC")) %>%
    # select key columns
    dplyr::select(SAMPLE, all_of(TARGET_PHENO)) %>%
    # ensure that the phenotype column is numeric
    dplyr::mutate(dplyr::across(all_of(TARGET_PHENO),
                                ~ as.numeric(.x))) %>%
    # randomise phenotype
    dplyr::mutate(dplyr::across(all_of(TARGET_PHENO),
                                ~ sample(.x)))

## Filter genotypes for those that have phenotypes
in_list[["genotypes"]] = in_list[["genotypes"]] %>%
    dplyr::filter(in_list[["sample_order"]] %in% phenos$SAMPLE)

## Filter and order phenotypes
in_list[["phenotypes"]] = phenos %>%
    # filter phenotypes for those with genotypes
    dplyr::filter(SAMPLE %in% in_list[["sample_order"]]) %>%
    # join to `sample_order` to ensure phenotypes are in the correct order   
    dplyr::left_join(tibble::tibble(SAMPLE = in_list[["sample_order"]]),
                     .,
                     by = "SAMPLE") %>%
    # remove NAs (created by the samples that have genotypes but not phenotypes)
    tidyr::drop_na() %>%
    # the GridLMM code doesn't work with tibbles
    as.data.frame()
            
# Run GWAS

out = run_gwas(d = in_list[["genotypes"]],
               m = in_list[["positions"]],
               p = in_list[["phenotypes"]]
              )

# Write results to file

out = saveRDS(out, OUT_FILE)

