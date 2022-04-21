# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(biomaRt)

# Get variables

## Debug
#RES = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/true/hdrr/None/5000/0.8/intercept.loco.mlma"
#MIN_P = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco_invnorm/min_p/hdrr/None/5000/0.8/intercept.csv"
#SITES = "/hps/nobackup/birney/users/ian/somites/sites_files/F0_Cab_Kaga/hdrr/homo_divergent/F1_het_min_DP.txt"
#BIN_LENGTH = "5000" %>% 
#  as.numeric()

## True
RES = snakemake@input[["res"]]
MIN_P = snakemake@input[["min_p"]]
SITES = snakemake@input[["sites"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT_BINS_ALL = snakemake@output[["bins_all"]]
OUT_BINS_UNQ = snakemake@output[["bins_unq"]]
OUT_SNPS = snakemake@output[["snps"]]
#OUT_SNPS_UNQ = snakemake@output[["snps_unq"]]

# Read in files

min_p = readr::read_csv(MIN_P) %>% 
  dplyr::pull(MIN_P) %>% 
  min(.)

results = readr::read_tsv(RES) %>% 
  # Add BIN_START and _END
  dplyr::mutate(BIN_START = (bp * BIN_LENGTH) + 1,
                BIN_END = (bp + 1) * BIN_LENGTH) %>% 
  # Pull significant loci
  dplyr::filter(p < min_p) %>% 
  # Separate SNP into CHROM and LOC
  tidyr::separate(SNP, into = c(NA, "BIN"), sep = ":") %>% 
  dplyr::mutate(BIN = as.numeric(BIN)) %>% 
  # Get columns in correct order
  dplyr::select(CHROM = Chr, BIN, BIN_START, BIN_END, everything(), -bp)


sites = readr::read_tsv(SITES,
                        col_names = c("CHROM", "POS_1", "POS_2", "REF", "ALT", "Cab_GT", "Kaga_GT")) %>% 
  # get Cab and Kaga alleles
  dplyr::mutate(CAB_ALLELE = dplyr::case_when(Cab_GT == "0/0" ~ REF,
                                              Cab_GT == "1/1" ~ ALT),
                KAGA_ALLELE = dplyr::case_when(Kaga_GT == "0/0" ~ REF,
                                               Kaga_GT == "1/1" ~ ALT)) %>% 
  dplyr::select(CHROM, POS = POS_1, REF, ALT, CAB_ALLELE, KAGA_ALLELE) %>% 
  # get bin
  dplyr::mutate(BIN = floor(POS / BIN_LENGTH),
                BIN_START = (BIN * BIN_LENGTH) + 1,
                BIN_END = (BIN + 1) * BIN_LENGTH) %>% 
  # bind with results to pull the SNPs within the significant loci
  dplyr::right_join(., results %>% 
                      dplyr::select(CHROM, BIN, A1, A2, Freq, b, se, p),
                    by = c("CHROM", "BIN"))

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


# BINS

## convert hits to genomic ranges
bin_loc_r = results %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM",
                                          start.field = "BIN_START",
                                          end.field = "BIN_END",
                                          ignore.strand = T)


## find overlaps
olaps_bins = GenomicRanges::findOverlaps(bin_loc_r, olat_genes_r)

# Pull out data frame of hits

bin_hits_all = dplyr::bind_cols(results %>% 
                              dplyr::slice(olaps_bins@from),
                              olat_genes[olaps_bins@to, ])

bin_hits_unq = olat_genes[unique(olaps_bins@to), ]


## SNPs
#
### convert hits to genomic ranges
#snp_loc_r = sites %>% 
#  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM",
#                                          start.field = "POS",
#                                          end.field = "POS",
#                                          ignore.strand = T)
#
#
### find overlaps
#olaps_snps = GenomicRanges::findOverlaps(snp_loc_r, olat_genes_r)
#
## Pull out data frame of hits
#
#snp_hits_all = dplyr::bind_cols(sites %>% 
#                                  dplyr::slice(olaps_snps@from),
#                                olat_genes[olaps_snps@to, ])
#
#snp_hits_unq = olat_genes[unique(olaps_snps@to), ]


# Write to files

readr::write_csv(bin_hits_all, OUT_BINS_ALL)
readr::write_csv(bin_hits_unq, OUT_BINS_UNQ)
readr::write_csv(sites, OUT_SNPS)
#readr::write_csv(snp_hits_unq, OUT_SNPS_UNQ)
