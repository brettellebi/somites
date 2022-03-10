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


```r
COUNTER = 0
lapply(results_list[1], function(RES){
  COUNTER <<- COUNTER + 1
  # Get bin length
  BIN_LENGTH = RES[["BIN_LENGTH"]] %>% 
    as.numeric()
  
  # Clean data frame
  test_results = RES[["RESULTS"]] %>% 
    dplyr::left_join(med_chr_lens, by = c("Chr" = "chr")) %>% 
    # add x-coord
    dplyr::mutate(X_COORD = pos + TOT) %>% 
    # change column names
    dplyr::rename(CHROM = Chr, BIN_START = pos) %>% 
    # add BIN_END
    dplyr::mutate(BIN_END = BIN_START + BIN_LENGTH - 1) %>% 
    # add locus
    dplyr::mutate(LOCUS = paste(CHROM, BIN_START, sep = ":")) %>%
    # target or not
    dplyr::mutate(TARGET = dplyr::if_else(LOCUS %in% sample_genos_list[[RES[["BIN_LENGTH"]]]]$LOCUS,
                                          "yes",
                                          "no"),
                  TARGET = factor(TARGET, levels = c("yes", "no"))) %>% 
    # create vector of colours
    dplyr::mutate(COLOUR = dplyr::case_when(TARGET == "yes" ~ names(gwas_pal)[1],
                                            gtools::even(CHROM) ~ names(gwas_pal)[2],
                                            gtools::odd(CHROM) ~ names(gwas_pal)[3]),
                  # order so that `target` is plotted last, at the front
                  COLOUR = factor(COLOUR, levels = rev(names(gwas_pal))),
                  SHAPE = dplyr::if_else(TARGET == "yes",
                                         18,
                                         20),
                  SIZE = dplyr::if_else(TARGET == "yes",
                                         1,
                                         0.5),
                  ALPHA = dplyr::if_else(TARGET == "yes",
                                         1,
                                         0.5)
                  )
  
  # Plot
  p1 = test_results %>% 
    ggplot(aes(x = X_COORD,
               y = -log10(p_value_REML),
               colour = COLOUR,
               shape = SHAPE,
               size = SIZE,
               alpha = ALPHA,
               label = CHROM,
               label2 = BIN_START,
               label3 = BIN_END)) + 
    geom_point() +
    aes(group = rev(TARGET)) +
    scale_color_manual(values = gwas_pal) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_alpha_identity() +
    scale_x_continuous(breaks = med_chr_lens$MID_TOT, 
                       labels = med_chr_lens$chr) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) +
    guides(colour = "none") + 
    ggtitle(paste("Bin length:", BIN_LENGTH)) +
    xlab("Chromosome") +
    ylab("-log10(p-value)")
  
  out = ggplotly(p1, tooltip = c("CHROM", "BIN_START", "BIN_END"))
  
  return(out)
})
#> $`20000`
```

