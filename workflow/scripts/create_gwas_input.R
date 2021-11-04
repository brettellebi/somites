# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

GENO_FILE = snakemake@input[["genotypes"]]
PHENO_FILE = snakemake@input[["phenotypes"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>%
    as.numeric()
OUT_FILE = snakemake@output[[1]]

# Read in and wrangle haplotype block data

df = readr::read_tsv(GENO_FILE,
                     col_types = "ciiidii") %>% 
    dplyr::mutate(SAMPLE = sample %>%
                    basename(.) %>% 
                    stringr::str_remove(".txt") %>% 
                    as.numeric(.),
                  BIN_START = (bin - 1) * BIN_LENGTH + 1,
                  BIN_END = bin * BIN_LENGTH) %>% 
    # recode state to make 0 == "Cab"
    dplyr::mutate(STATE = dplyr::recode(state,
                                        `0` = 2,
                                        `1` = 1,
                                        `2` = 0)) %>% 
    dplyr::select(SAMPLE, CHROM = chr, BIN = bin, BIN_START, BIN_END, STATE) %>%
    # order by SAMPLE
    dplyr::arrange(SAMPLE, CHROM, BIN)

# Filter for loci with > 1 genotype across all samples

## Widen data frame
gt_df = df %>%
    tidyr::pivot_wider(names_from = SAMPLE, values_from = STATE)

## Pull out sample names
sample_names = colnames(gt_df)[5:ncol(gt_df)]

## Pull out matrix of genotypes
gt_mat = as.matrix(gt_df[, 5:ncol(gt_df)])

## Get indexes of loci with > 1 genotype
bins_to_keep = logical()

for (ROW in 1:nrow(gt_mat)){
    ## get unique values in each row
    out = unique(gt_mat[ROW, ])
    ## remove NAs
    out = out[!is.na(out)]
    # if more than one value, return TRUE
    if (length(out) > 1) {
        bins_to_keep[ROW] = TRUE
    }
    # if just one value (i.e. if all samples are the same genotype at that locus), return false 
    else {
        bins_to_keep[ROW] = FALSE
    }
}

## filter gt_final
gt_filt = gt_df %>% 
    dplyr::filter(bins_to_keep) %>% 
    # recode genotypes to -1, 0, 1
    dplyr::mutate(dplyr::across(-c("CHROM", "BIN", "BIN_START", "BIN_END"),
                                ~dplyr::recode(.x,
                                               `0` = -1,
                                               `1` = 0,
                                               `2` = 1))) %>% 
    # order
    dplyr::arrange(CHROM, BIN_START)  

# Create final input for GWAS as list

out_list = list()

## Genotypes
out_list[["genotypes"]] = gt_filt %>% 
    dplyr::select(-c(CHROM, BIN, BIN_START, BIN_END)) %>% 
    # convert to matrix
    as.matrix(.) %>% 
    # transpose to put samples as rows
    t(.) %>% 
    # convert to data frame (GridLMM code doesn't work with tibbles)
    as.data.frame(.)

## Sample order
out_list[["sample_order"]] = sample_names

## Positions
out_list[["positions"]] = gt_filt %>% 
    dplyr::select(CHROM, BIN_START, BIN_END) %>%
    # The GridLMM code doesn't work with tibbles
    as.data.frame(.)

# Save to file

saveRDS(out_list, OUT_FILE)
