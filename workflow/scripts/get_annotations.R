# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(here)
library(tidyverse)
library(biomaRt)

# Get variables

## Debug
#DATE_OF_ASSOC_TEST = "20220214"
#SITE_FILTER = "all_sites"
#COVARIATES = "Microscope"
#INVERSE_NORM = "TRUE"
#BIN_LENGTH = 5000
#OUT_DIR = here::here("data", DATE_OF_ASSOC_TEST, "annotations")
#TARGET_PHENO = "unsegmented_psm_area"
#
#IN_FILE = file.path("/hps/nobackup/birney/users/ian/somites/association_testing",
#                    DATE_OF_ASSOC_TEST,
#                    FILTER,
#                    "true_results",
#                    TARGET_PHENO,
#                    COVARIATES,
#                    INVERSE_NORM,
#                    paste(BIN_LENGTH, ".rds", sep = ""))
#PERMS_FILES = list.files(file.path("/hps/nobackup/birney/users/ian/somites/association_testing",
#                                   DATE_OF_ASSOC_TEST,
#                                   FILTER,
#                                   "permutations",
#                                   TARGET_PHENO,
#                                   COVARIATES,
#                                   INVERSE_NORM,
#                                   BIN_LENGTH),
#                         full.names = T)


## True
GWAS_RESULTS = snakemake@input[["gwas_results"]]
PERMS_FILES = snakemake@input[["perm_results"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT_FILE = snakemake@output[["csv"]]


# Read in files

## Results
RESULTS = readRDS(GWAS_RESULTS)

if ("p_value_REML" %in% colnames(RESULTS$results)){
  P_COL = "p_value_REML"
} else {
  P_COL = "p_value_REML.1"
}
### Rename column in results
if (P_COL == "p_value_REML.1"){
  RESULTS$results = RESULTS$results %>% 
    dplyr::rename(p_value_REML = p_value_REML.1)
}

## Permutations
PERM_LIST = purrr::map(PERMS_FILES, readRDS)

perm_df = purrr::map_dfr(PERM_LIST, function(PERM){
  OUT = tibble::tibble(MIN_P = PERM$results %>% 
                         dplyr::select(dplyr::all_of(P_COL)) %>%
                         min(., na.rm = T)
  )
}, .id = "SEED")

# Get minimum
SIG_LEVEL = min(perm_df$MIN_P)

# Pull significant loci

SIG_LOCS = RESULTS$results %>% 
  dplyr::filter(p_value_REML < SIG_LEVEL) %>% 
  dplyr::select(CHROM = Chr,
                BIN_START = pos) %>% 
  dplyr::mutate(BIN_END = BIN_START + BIN_LENGTH - 1)

# Get annotations

## Select dataset
olat_mart = biomaRt::useEnsembl(biomart = "ensembl", dataset = "olatipes_gene_ensembl", mirror = "uswest")

olat_attr = biomaRt::listAttributes(olat_mart)

olat_genes = biomaRt::getBM(attributes = c("chromosome_name",
                                           "start_position",
                                           "end_position",
                                           "ensembl_gene_id",
                                           "hgnc_symbol",
                                           "ensembl_exon_id",
                                           "description",
                                           "strand",
                                           "transcript_start",
                                           "transcript_end"),
                            mart = olat_mart) 

olat_genes_r = olat_genes %>% 
  # change strand characters
  dplyr::mutate(strand = dplyr::recode(.$strand,
                                       `1` = "+",
                                       `-1` = "-")
  ) %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "chromosome_name",
                                          start.field = "start_position",
                                          end.field = "end_position")

## convert hits to genomic ranges
sig_loc_r = SIG_LOCS %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM",
                                          start.field = "BIN_START",
                                          end.field = "BIN_END",
                                          ignore.strand = T)


## find overlaps
olaps = GenomicRanges::findOverlaps(sig_loc_r, olat_genes_r)

# Pull out data frame of hits

hits = olat_genes[unique(olaps@to), ]

# Save to file

readr::write_csv(hits,
                 OUT_FILE)
