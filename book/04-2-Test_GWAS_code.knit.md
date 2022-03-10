---
zotero: PhD
---

# Association testing with simulated phenotypes

The code for simulating phenotypes and using them to test the GWLS code is set out in this Snakefile, specifically rules `simulate_phenotypes` and `test_gwls`: https://github.com/brettellebi/somites/blob/master/workflow/rules/07_assocation_testing.smk


```r
date_of_assoc_test = "20220214"
```



```r
library(here)
source(here::here("book/source/04-Association_testing.R"))
```

## Read in F2 recombination blocks

### Read data


```r
in_dir = "/nfs/research/birney/users/ian/somites/recombination_blocks/F2/all_sites"

in_files = list.files(in_dir, full.names = T)

# Read into list
data_list = purrr::map(in_files, function(FILE){
  out = readr::read_tsv(FILE,
                        col_types = "ciiidii")
})
# Set names as bin length
names(data_list) = basename(in_files) %>% 
  stringr::str_remove(".txt")
# Reorder
data_list = data_list[order(as.numeric(names(data_list)))]

counter = 0
df_list = purrr::map(data_list, function(data){
  counter <<- counter + 1
  # set bin length
  bin_length = as.numeric(names(data_list)[counter])
  # add bin start and end coordinates
  df = data %>% 
    dplyr::mutate(SAMPLE = basename(sample) %>% 
                    stringr::str_remove(".txt") %>% 
                    as.numeric(.),
                  BIN_START = (bin - 1) * bin_length + 1,
                  BIN_END = bin * bin_length) %>% 
    # recode state to make 0 == "Cab"
    dplyr::mutate(STATE = dplyr::recode(state,
                                        `0` = 2,
                                        `1` = 1,
                                        `2` = 0)) %>% 
    dplyr::select(SAMPLE, CHROM = chr, BIN = bin, BIN_START, BIN_END, STATE)
  
  return(df)
})

```

### How many possible blocks have at least one call? 

#### Count total existing bins


```r
bin_lengths = as.integer(names(df_list))
names(bin_lengths) = bin_lengths

n_bins = purrr::map(bin_lengths, function(BIN_LENGTH){
  purrr::map_int(med_chr_lens$end, function(CHR_END){
    out = seq(from = 1, to = CHR_END, by = BIN_LENGTH)
    
    return(length(out))
  })
})

n_bins_total = purrr::map_int(n_bins, sum)
```

#### Proportion of total bins with calls


```r
# Get total number of samples
n_samples = df_list$`5000`$SAMPLE %>%
  unique(.) %>% 
  length(.)

# Build DF of bins
n_bins_df = purrr::map_dfr(1:length(df_list), function(COUNTER){
  # Bin length
  bin_length = as.numeric(names(df_list)[COUNTER])
  # Number of total bins
  n_bins = n_bins_total[COUNTER]
  # Number of bins with calls
  n_bins_with_calls = df_list[[COUNTER]] %>% 
    dplyr::distinct(CHROM, BIN_START) %>% 
    nrow(.)
  # Number of bins with calls for each 
  n_bins_no_missing = df_list[[COUNTER]] %>% 
    dplyr::count(CHROM, BIN_START) %>%
    dplyr::filter(n == n_samples) %>%
    nrow(.)
  
  # Build final data frame
  out = tibble::tibble(BIN_LENGTH = bin_length,
                       N_BINS = n_bins,
                       N_BINS_WITH_CALLS = n_bins_with_calls,
                       N_BINS_NO_MISSING = n_bins_no_missing) %>% 
    dplyr::mutate(PROP_BINS_WITH_CALLS = N_BINS_WITH_CALLS / N_BINS,
                  PROP_BINS_NO_MISSING = N_BINS_NO_MISSING / N_BINS )
  
  return(out)
})

DT::datatable(n_bins_df)
```

```{=html}
<div id="htmlwidget-fd228a3db6c0125fafe9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fd228a3db6c0125fafe9">{"x":{"filter":"none","data":[["1","2","3","4"],[5000,10000,15000,20000],[146819,73414,48947,36712],[111760,62718,43880,33843],[19,50,81,142],[0.761209380257324,0.854305718255374,0.896479865977486,0.921851165831336],[0.000129411043529788,0.000681069005911679,0.00165485116554641,0.0038679450860754]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>BIN_LENGTH<\/th>\n      <th>N_BINS<\/th>\n      <th>N_BINS_WITH_CALLS<\/th>\n      <th>N_BINS_NO_MISSING<\/th>\n      <th>PROP_BINS_WITH_CALLS<\/th>\n      <th>PROP_BINS_NO_MISSING<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

### Merge recombination blocks

#### List with final recoded genotypes (including NAs)


```r
gt_list = purrr::map(df_list, function(BIN_LENGTH_DF){
  
  # Widen data frame
  gt_final = BIN_LENGTH_DF %>% 
    tidyr::pivot_wider(names_from = SAMPLE, values_from = STATE)
  
  # Pull out matrix of genotypes
  gt_mat = as.matrix(gt_final[, 5:ncol(gt_final)])
  
  # Get indexes of loci with > 1 genotype
  bins_to_keep = logical()
  
  for (ROW in 1:nrow(gt_mat)){
    # get unique values in each row
    out = unique(gt_mat[ROW, ])
    # remove NAs
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
  
  # filter gt_final
  gt_filt = gt_final %>% 
    dplyr::filter(bins_to_keep) %>% 
    # recode genotypes to -1, 0, 1
    dplyr::mutate(dplyr::across(-c("CHROM", "BIN", "BIN_START", "BIN_END"),
                                ~dplyr::recode(.x,
                                               `0` = -1,
                                               `1` = 0,
                                               `2` = 1))) %>% 
    # order
    dplyr::arrange(CHROM, BIN_START)  
  
  return(gt_filt)
}) 

# Show first 10 columns
gt_list$`20000`[, 1:10] %>% 
  head(.) %>% 
  DT::datatable(.) 
```

```{=html}
<div id="htmlwidget-22cf139c18731560b797" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-22cf139c18731560b797">{"x":{"filter":"none","data":[["1","2","3","4","5","6"],[1,1,1,1,1,1],[1,2,3,4,5,6],[1,20001,40001,60001,80001,100001],[20000,40000,60000,80000,100000,120000],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>CHROM<\/th>\n      <th>BIN<\/th>\n      <th>BIN_START<\/th>\n      <th>BIN_END<\/th>\n      <th>1<\/th>\n      <th>2<\/th>\n      <th>4<\/th>\n      <th>6<\/th>\n      <th>8<\/th>\n      <th>9<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

#### List with final recoded genotypes (complete cases only)


```r
gt_nomiss_list = purrr::map(gt_list, function(BIN_LENGTH_DF){
  BIN_LENGTH_DF %>%
    dplyr::filter(complete.cases(.))
})
```


## Simulate phenotype

### Extract samples of genotypes


```r
# Set directory
out_dir_test = file.path(lts_dir, "association_testing/20211027_test")
```

These must be written to a file because `PhenotypeSimulator` reads delimited genotypes from files.


```r
# Get random 10 loci
set.seed(10)
n_loci = 10
# NOTE: PhenotypeSimulator::readStandardGenotypes states that the genotype file must
# delim: a [delimter]-delimited file of [(NrSNPs+1) x (NrSamples+1)] genotypes with the snpIDs in the first column and the sampleIDs in the first row and genotypes encoded as numbers between 0 and 2 representing the (posterior) mean genotype, or dosage of the minor allele.

sample_gts = purrr::map(seq_along(gt_nomiss_list), function(COUNTER){
  
  # Pull out random SNPs
  snp_sample = gt_nomiss_list[[COUNTER]] %>% 
    dplyr::slice_sample(n = n_loci) %>% 
    dplyr::arrange(CHROM, BIN_START) %>% 
    # create locus column
    dplyr::mutate(LOCUS = paste(CHROM, BIN_START, sep = ":")) %>% 
    # remove superfluous columns
    dplyr::select(-c(CHROM, BIN, BIN_START, BIN_END)) %>% 
    # recode back to 0,1,2
    dplyr::mutate(dplyr::across(-LOCUS,
                                ~dplyr::recode(.x,
                                               `-1` = 0,
                                               `0` = 1,
                                               `1` = 2))) %>% 
    # reorder columns
    dplyr::select(LOCUS, everything())
  
})
names(sample_gts) = names(gt_nomiss_list)

purrr::map(seq_along(sample_gts), function(COUNTER){
  # save to file
  bin_length = names(sample_gts)[COUNTER]
  out_file = file.path(out_dir_test, "target_loci", paste(bin_length, ".csv", sep = ""))
  readr::write_csv(sample_gts[[COUNTER]], out_file)
})

saveRDS(sample_gts, file.path(out_dir_test, "target_loci/sample_gts.rds"))
```




### Simulate phenotype


```r
sim_path = file.path(out_dir_test, "simulated_phenotypes/20211103_sim_phenos.rds")
```


```r
set.seed(5671)
# N samples
N = n_samples
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

sim_phenos = purrr::map(seq_along(sample_gts), function(COUNTER){
  # get sample file path
  bin_length = names(sample_gts)[COUNTER]
  gt_sample_file = file.path(out_dir_test, "target_loci", paste(bin_length, ".csv", sep = ""))
    
  sim_pheno = PhenotypeSimulator::runSimulation(N = N, P = P, 
                                                genVar = genVar, h2s = h2s, 
                                                noiseVar = noiseVar, phi = phi,
                                                cNrSNP = n_loci,
                                                genotypefile = gt_sample_file,
                                                format = "delim",
                                                genoDelimiter = ",")
  
  return(sim_pheno)
})
names(sim_phenos) = names(sample_gts)

saveRDS(sim_phenos, sim_path)

# Write as .xlsx to use in same Snakemake code as true GWLS
lapply(seq_along(sim_phenos), function(COUNTER){
  out = tibble::tibble(fish = colnames(sample_gts[[COUNTER]])[-1],
                       Y = sim_phenos[[COUNTER]]$phenoComponentsFinal$Y)
  # set path for output
  out_file = file.path(lts_dir, 
                       "association_testing/20211027_test/simulated_phenotypes",
                       paste(names(sim_phenos)[COUNTER], ".xlsx", sep = ""))
  # write to file
  openxlsx::write.xlsx(out, out_file, overwrite = T)
})
```



## Reformat genotypes for GridLMM

### Complete cases


```r
final_nomiss = purrr::map(seq_along(gt_nomiss_list), function(COUNTER){
  out = list()
  
  # Genotypes
  out[["genotypes"]] = gt_nomiss_list[[COUNTER]] %>% 
    dplyr::select(-c(CHROM, BIN, BIN_START, BIN_END)) %>% 
    # convert to matrix
    as.matrix(.) %>% 
    # transpose to put samples as rows
    t(.) %>% 
    # convert to data frame
    as.data.frame(.)
  
  # Positions
  out[["positions"]] = gt_nomiss_list[[COUNTER]] %>% 
    dplyr::select(CHROM, BIN_START, BIN_END)
  
  # Phenotypes
  out[["phenotypes"]] = data.frame(fish = rownames(sim_phenos[[COUNTER]]$phenoComponentsFinal$Y),
                                   Y = sim_phenos[[COUNTER]]$phenoComponentsFinal$Y %>%
                                     as.numeric()
                                   )
    
  return(out)
})
names(final_nomiss) = names(gt_nomiss_list)
```

### With NAs


```r
final_wimiss = purrr::map(seq_along(gt_list), function(COUNTER){
  out = list()
  
  # Genotypes
  out[["genotypes"]] = gt_list[[COUNTER]] %>% 
    dplyr::select(-c(CHROM, BIN, BIN_START, BIN_END)) %>% 
    # convert to matrix
    as.matrix(.) %>% 
    # transpose to put samples as rows
    t(.) %>% 
    # convert to data frame
    as.data.frame(.)
  
  # Positions
  out[["positions"]] = gt_list[[COUNTER]] %>% 
    dplyr::select(CHROM, BIN_START, BIN_END)
  
  # Phenotypes
  out[["phenotypes"]] = data.frame(SAMPLE = rownames(sim_phenos[[COUNTER]]$phenoComponentsFinal$Y),
                                   Y = sim_phenos[[COUNTER]]$phenoComponentsFinal$Y %>%
                                     as.numeric()
                                   )
    
  return(out)
})
names(final_wimiss) = names(gt_list)
```

## Read in sampled genotypes


```r
SAMPLE_GTS_DIR = "/nfs/research/birney/users/ian/somites/association_testing/20220118"

# List files
TARGET_DIRS = list.files(SAMPLE_GTS_DIR, pattern = "sample_genos", full.names = T, recursive = T, include.dirs = T)
TARGET_FILES = list.files(TARGET_DIRS, full.names = T, recursive = T); names(TARGET_FILES) = TARGET_FILES

sample_genos_list = purrr::map_dfr(TARGET_FILES, function(FILE){
  readr::read_csv(FILE)
}, .id = "FILENAME")
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 570
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (569): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 10 Columns: 548
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr   (1): LOCUS
#> dbl (547): 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


## Read in GWLS results


```r
RESULTS_DIR = "/nfs/research/birney/users/ian/somites/association_testing/20220214"

# List files
TARGET_DIRS = list.files(RESULTS_DIR, pattern = "test_results", full.names = T, recursive = T, include.dirs = T)
TARGET_FILES = list.files(TARGET_DIRS, full.names = T, recursive = T); names(TARGET_FILES) = TARGET_FILES

# Create DF to separate file name into variables
meta_df = tibble::tibble(FILENAME = TARGET_FILES) %>% 
  tidyr::separate(col = FILENAME,
                  sep = "/",
                  into = c(rep(NA, 9), "SITE_FILTER", NA, "BIN_LENGTH")) %>% 
  # remove extension from BIN_LENGTH
  dplyr::mutate(BIN_LENGTH = stringr::str_remove(BIN_LENGTH, ".rds"))

# Get indexes of target combinations
TARGET_INDEXES = which(meta_df$SITE_FILTER %in% c("all_sites", "no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop") & meta_df$BIN_LENGTH %in% c(5000, 20000))

# Read into list
results_list = purrr::map(seq_along(TARGET_FILES), function(INDEX){
  OUT = list()
  
  # Add metadata info
  OUT[["SITE_FILTER"]] = meta_df %>% 
    dplyr::slice(INDEX) %>% 
    dplyr::pull("SITE_FILTER")
  OUT[["BIN_LENGTH"]] = meta_df %>% 
    dplyr::slice(INDEX) %>% 
    dplyr::pull("BIN_LENGTH")
  
  # Read in results
  results = readRDS(TARGET_FILES[INDEX])
  ## Add to list
  OUT[["RESULTS"]] = results$results
  
  return(OUT)
})

# Filter for desired extreme combinations
list_to_plot = results_list[TARGET_INDEXES]
```

## Plot





























