---
zotero: PhD
---

# Association testing with simulated phenotypes

The code for simulating phenotypes and using them to test the GWLS code is set out in this Snakefile, specifically rules `simulate_phenotypes` and `test_gwls`: https://github.com/brettellebi/somites/blob/master/workflow/rules/07_association_testing.smk


```r
date_of_assoc_test = "20220214"
RECOM_BLOCKS_DIR = "/nfs/research/birney/users/ian/somites/recombination_blocks/F2/all_sites"
SAMPLE_GTS_DIR = "/hps/nobackup/birney/users/ian/somites/association_testing/20220214/all_sites/sample_genos"
RESULTS_DIR = "/hps/nobackup/birney/users/ian/somites/association_testing/20220214"
```



```r
library(here)
source(here::here("book/source/04-Association_testing.R"))
```

## Read in sampled genotypes


```r
# List files
TARGET_FILES = list.files(SAMPLE_GTS_DIR, full.names = T) ; names(TARGET_FILES) = TARGET_FILES

sample_genos_list = purrr::map(TARGET_FILES, function(FILE){
  readr::read_csv(FILE) %>% 
    dplyr::select(LOCUS)
})
names(sample_genos_list) = TARGET_FILES %>% 
  basename() %>% 
  stringr::str_remove(".csv")
```


## Read in GWLS results


```r
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

names(results_list) = TARGET_FILES %>% 
  basename() %>% 
  stringr::str_remove(".rds")
```

## Plot



