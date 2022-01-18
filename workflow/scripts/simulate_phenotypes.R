# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(data.table)
library(PhenotypeSimulator)
library(openxlsx)

# Get variables

## Debugging
IN_FILE = "/nfs/research/birney/users/ian/somites/association_testing/20220118_true/no_repeat_sites/inputs/5000.rds"
OUT_SAMPLE_GENOS = "/nfs/research/birney/users/ian/somites/association_testing/20220118_test/no_repeat_sites/sample_genos/5000.csv"
OUT_SIM_PHENOS = "/nfs/research/birney/users/ian/somites/association_testing/20220118_test/no_repeat_sites/sim_phenos/5000.xlsx"
N_SAMPLE_GTS = 10

## True
IN_FILE = snakemake@input[[1]]
OUT_SAMPLE_GENOS = snakemake@output[["sample_genos"]]
OUT_SIM_PHENOS = snakemake@output[["sim_phenos"]]
N_SAMPLE_GTS = snakemake@params[["n_sample_gts"]]

# Read in genotypes file and create sample GTs file

### NOTE: These must be written to a file because `PhenotypeSimulator` reads delimited genotypes from files.

# NOTE: PhenotypeSimulator::readStandardGenotypes states that the genotype file must
# delim: a [delimter]-delimited file of [(NrSNPs+1) x (NrSamples+1)] genotypes with the snpIDs in the first column and the sampleIDs in the first row and genotypes encoded as numbers between 0 and 2 representing the (posterior) mean genotype, or dosage of the minor allele.

## Read in genotypes file
input = readRDS(IN_FILE)
## Create data frame in format required for PhenotypeSimulator
in_df = cbind(
  # Add position
  paste(input$positions$CHROM, input$positions$BIN_START, sep = ":"),
  # Add genotypes
  data.table::transpose(input$genotypes)
  )
## Add sample names
colnames(in_df) = c("LOCUS", input$sample_order)
## Get random 10 loci
set.seed(10)
## Select sample and write to file
in_df = in_df %>% 
  # Get complete cases
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::slice_sample(n = N_SAMPLE_GTS) %>% 
  ## recode back to 0,1,2
  dplyr::mutate(dplyr::across(-LOCUS,
                              ~dplyr::recode(.x,
                                             `-1` = 0,
                                             `0` = 1,
                                             `1` = 2)))
## Write to file
in_df %>% 
  readr::write_csv(OUT_SAMPLE_GENOS)

# Simulate phenotype

set.seed(5671)
# Change number of target SNPs if there are not enough SNPs with complete cases
if (N_SAMPLE_GTS > nrow(in_df)){
  N_SAMPLE_GTS = nrow(in_df)
}
# N samples
N = length(input$sample_order)
# N phenotypes
P = 1
# Proportion of total genetic variance
genVar = 0.5
# Proportion of genetic variance of genetic variant effects
h2s = 1
# Proportion of total noise variance
noiseVar = 0.5
# Proportion of noise variance of observational noise effects
phi = 1

sim_pheno = PhenotypeSimulator::runSimulation(N = N, P = P, 
                                              genVar = genVar, h2s = h2s, 
                                              noiseVar = noiseVar, phi = phi,
                                              cNrSNP = N_SAMPLE_GTS,
                                              genotypefile = OUT_SAMPLE_GENOS,
                                              format = "delim",
                                              genoDelimiter = ",")
  

# Write as .xlsx to use in same Snakemake code as true GWLS
out = tibble::tibble(fish = input$sample_order,
                     Y = sim_pheno$phenoComponentsFinal$Y)

# write to file
openxlsx::write.xlsx(out, OUT_SIM_PHENOS, overwrite = T)
