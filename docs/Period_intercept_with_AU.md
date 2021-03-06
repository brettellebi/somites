# Association testing with only AU microscope

## Load libraries and variables


```r
library(tidyverse)
library(GridLMM)
library(KRLS)

GENO_FILE = "/hps/nobackup/birney/users/ian/somites/association_testing/20220321/all_sites/inputs/5000.rds"
PHENO_FILE = here::here("data/20220321_phenotypes.xlsx") # True phenotypes
GWLS_SOURCE_FILE = here::here("workflow/scripts/run_gwas_source.R")
MANHAT_SOURCE_FILE = here::here("workflow/scripts/get_manhattan_source.R")
BIN_LENGTH = 5000
TARGET_PHENO = "intercept"
MICROSCOPE = "AU"
PLOT_DIR = here::here("book/plots/20220321/microscope_test")
ALPHA = 0.05

# Get GWAS functions

source(GWLS_SOURCE_FILE)
source(MANHAT_SOURCE_FILE)
```

## Load genotypes and positions


```r
in_list = readRDS(GENO_FILE)
```

## Read in phenotypes

```r


## Read in file and wrangle
phenos = readxl::read_xlsx(PHENO_FILE) %>%
    # adjust sample names
    dplyr::mutate(SAMPLE = fish %>% stringr::str_remove("KC")) %>%
    # select key columns
    dplyr::select(SAMPLE, all_of(TARGET_PHENO), Microscope) %>%
    # ensure that the phenotype column is numeric
    dplyr::mutate(dplyr::across(all_of(TARGET_PHENO),
                                ~ as.numeric(.x))) %>% 
    # THIS IS THE (big) CHANGE: filter for microscope
    dplyr::filter(Microscope == MICROSCOPE)
```

## Filter for samples with both genos and phenos


```r
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

## Filter genotypes for those that have phenotypes
in_list[["genotypes"]] = in_list[["genotypes"]] %>%
    dplyr::slice(in_list[["sample_order"]] %in% in_list[["phenotypes"]]$SAMPLE %>% 
                   which())

## Filter sample_order for those that have phenotypes
in_list[["sample_order"]] = in_list[["phenotypes"]]$SAMPLE

## Get number of samples
N_SAMPLES = in_list[["sample_order"]] %>% 
  length()
```

## Run GWAS


```r
# Run GWAS

out = run_gwas(d = in_list[["genotypes"]],
               m = in_list[["positions"]],
               p = in_list[["phenotypes"]],
               invers_norm = T
              )

saveRDS(out, paste("/hps/nobackup/birney/users/ian/somites/microscope_test/gwas_results/", MICROSCOPE, ".rds", sep = ""))
```


```r
out = readRDS(paste("/hps/nobackup/birney/users/ian/somites/microscope_test/gwas_results/", MICROSCOPE, ".rds", sep = ""))
```

## Run permutations

### Permute phenos


```r
seeds = 1:10

counter = 0
perm_phenos = purrr::map(seeds, function(SEED){
  counter <<- counter + 1
  # set seed
  set.seed(seeds[counter])
  # get original phenotypes
  orig_phenos = phenos
  # randomise
  phenos = orig_phenos %>% 
      # randomise phenotype
      dplyr::mutate(dplyr::across(all_of(TARGET_PHENO),
                                  ~ sample(.x)))
})
```

### Run GWLS


```r
perm_out = purrr::map(perm_phenos, function(PERM_PHENO){
  ## Get phenotypes
  phenos = PERM_PHENO
  
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
  
  ## Filter genotypes for those that have phenotypes
  in_list[["genotypes"]] = in_list[["genotypes"]] %>%
      dplyr::filter(in_list[["sample_order"]] %in% in_list[["phenotypes"]]$SAMPLE)
  
  ## Filter sample_order for those that have phenotypes
  in_list[["sample_order"]] = in_list[["phenotypes"]]$SAMPLE
              
  # Run GWAS
  
  out = run_gwas(d = in_list[["genotypes"]],
                 m = in_list[["positions"]],
                 p = in_list[["phenotypes"]],
                 invers_norm = T
                )
  
  return(out)
})

names(perm_out) = seeds

saveRDS(perm_out, paste("/hps/nobackup/birney/users/ian/somites/microscope_test/perms/", MICROSCOPE, ".rds", sep = ""))
```


```r
perm_out = readRDS(paste("/hps/nobackup/birney/users/ian/somites/microscope_test/perms/", MICROSCOPE, ".rds", sep = ""))
```

### Get minimum


```r
perm_df = purrr::map_dfr(perm_out, function(PERM){
  OUT = tibble::tibble(MIN_P = PERM$results$p_value_REML %>%
                         min(., na.rm = T)
  )
}, .id = "SEED")

# Get minimum
SIG_LEVEL = min(perm_df$MIN_P)

# Get bonferroni correction
SIG_BONF = ALPHA / ncol(in_list[["genotypes"]])
```

## Generate Manhattan plot


```r
out_clean = clean_gwas_res(out,
                           bin_length = BIN_LENGTH,
                           chr_lens = med_chr_lens)

# Plot
out_plot = plot_man(out_clean,
                    site_filter = "all_sites",
                    phenotype = TARGET_PHENO,
                    bin_length = BIN_LENGTH, 
                    gwas_pal = intercept_pal,
                    med_chr_lens = med_chr_lens,
                    sig_level = SIG_LEVEL,
                    bonferroni = SIG_BONF) +
                 labs(subtitle = paste("Microscope: ", MICROSCOPE, "\nCovariates: None\nn samples: ", N_SAMPLES, sep = ""))

out_plot
#> Warning: Removed 7 rows containing missing values
#> (geom_point).
```

<img src="Period_intercept_with_AU_files/figure-html/unnamed-chunk-11-1.png" width="672" />


```r
ggsave(file.path(PLOT_DIR, paste(MICROSCOPE, "_manhattan.png", sep = "")),
       out_plot,
       device = "png",
       width = 9.6,
       height = 6,
       units = "in",
       dpi = 400)
#> Warning: Removed 7 rows containing missing values
#> (geom_point).
```

