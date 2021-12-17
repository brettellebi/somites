---
zotero: PhD
---

# GWLS results


```r
library(here)
source(here::here("book/source/04-Association_testing.R"))
```

## Notes

* 20211104 association test performed on full sites file. Results here: `/nfs/research/birney/users/ian/somites/association_testing/20211104_true/results`
* 20211109 association test performed on sites excluding those overlapping repeat regions. Results here: `/nfs/research/birney/users/ian/somites/association_testing/20211109_true/results`

## Phenotype data

First 400 phenotype data: <https://github.com/brettellebi/somites/tree/master/data/20210917_First400_F2_DF.xlsx>.

>I also attach here the DataFrame of the phenotyping of the first 400 F2s (First400_F2_DF.xlsx)...

>Of interest for us for the association testing are the:
 intercept period -> intercept in the table
 mean period -> mean in the table
 
## Snakemake rules

Snakemake rules for running the GWAS over these phenotypes: <https://github.com/brettellebi/somites/blob/master/workflow/rules/07_assocation_testing.smk>

## Results

### All sites

#### Read in files


```r
target_phenos = c("mean", "intercept")
names(target_phenos) = target_phenos
bin_lengths = c(5000, 10000, 15000, 20000)
names(bin_lengths) = bin_lengths
results_dir = "/nfs/research/birney/users/ian/somites/association_testing/20211109_true/all_sites/results/"

gwas_true = purrr::map(target_phenos, function(PHENO){
  purrr::map(bin_lengths, function(BIN_LENGTH){
    readRDS(file.path(results_dir, PHENO, paste(BIN_LENGTH, ".rds", sep = "")))
  })
})
```

#### Manhattans

Process data

```r
plot_dat = purrr::map(seq_along(gwas_true), function(COUNTER_1){
  bin_list = purrr::map(seq_along(gwas_true[[COUNTER_1]]), function(COUNTER_2){
    BIN_LENGTH = names(gwas_true[[COUNTER_1]])[COUNTER_2] %>% 
      as.numeric()
    out = clean_gwas_res(gwas_true[[COUNTER_1]][[COUNTER_2]],
                         bin_length = BIN_LENGTH,
                         chr_lens = med_chr_lens)
    
    return(out)
  })
  names(bin_list) = names(gwas_true[[COUNTER_1]])
  
  return(bin_list)
})
names(plot_dat) = names(gwas_true)
```

Read in significance levels from permutations


```r
perm_file = here::here("data/20211109_permutation_mins.csv")
perms_df_mins = readr::read_csv(perm_file,
                                col_types = c("ccid"))
```

Plot

```r
plot_dir = here::here("book/plots/20211109/gwls_results/all_sites")
# Plot
lapply(seq_along(plot_dat), function(COUNTER_PHENO){
  lapply(seq_along(plot_dat[[COUNTER_PHENO]]), function(COUNTER_BINL){
    # Get pheno
    PHENO = names(plot_dat)[COUNTER_PHENO]
    # Get bin length
    BINL = names(plot_dat[[COUNTER_PHENO]])[COUNTER_BINL] %>% 
      as.numeric()
    # Get significant level
    SIG_LEVEL = perms_df_mins %>%
      dplyr::filter(SITE_FILTER == "all_sites", TARGET_PHENO == PHENO, BIN_LENGTH == BINL) %>% 
      dplyr::pull(MIN_P) %>% 
      -log10(.)
    
    # Plot
    out = plot_man(plot_dat[[COUNTER_PHENO]][[COUNTER_BINL]],
             phenotype = PHENO,
             bin_length = BINL,
             gwas_pal = gwas_pal[2:3],
             med_chr_lens = med_chr_lens,
             sig_line = SIG_LEVEL) +
      ylim(0,7)
    
    ggsave(file.path(plot_dir, paste(PHENO, "_", BINL, ".png", sep = "")),
           out,
           device = "png",
           width = 9.6,
           height = 6,
           units = "in",
           dpi = 400)
    
    out
  })
})
#> [[1]]
#> [[1]][[1]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```
#> 
#> [[1]][[2]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-2.png" width="672" />

```
#> 
#> [[1]][[3]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-3.png" width="672" />

```
#> 
#> [[1]][[4]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-4.png" width="672" />

```
#> 
#> 
#> [[2]]
#> [[2]][[1]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-5.png" width="672" />

```
#> 
#> [[2]][[2]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-6.png" width="672" />

```
#> 
#> [[2]][[3]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-7.png" width="672" />

```
#> 
#> [[2]][[4]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-5-8.png" width="672" />


### Excluding reads overlapping repeats

#### Read in files


```r
filter_type = "no_repeat_reads"
target_phenos = c("mean", "intercept")
names(target_phenos) = target_phenos
bin_lengths = c(5000, 10000, 15000, 20000)
names(bin_lengths) = bin_lengths
results_dir = file.path("/nfs/research/birney/users/ian/somites/association_testing/20211109_true", filter_type, "results")

gwas_true = purrr::map(target_phenos, function(PHENO){
  purrr::map(bin_lengths, function(BIN_LENGTH){
    readRDS(file.path(results_dir, PHENO, paste(BIN_LENGTH, ".rds", sep = "")))
  })
})
```

#### Manhattans

Process data

```r
plot_dat = purrr::map(seq_along(gwas_true), function(COUNTER_1){
  bin_list = purrr::map(seq_along(gwas_true[[COUNTER_1]]), function(COUNTER_2){
    BIN_LENGTH = names(gwas_true[[COUNTER_1]])[COUNTER_2] %>% 
      as.numeric()
    out = clean_gwas_res(gwas_true[[COUNTER_1]][[COUNTER_2]],
                         bin_length = BIN_LENGTH,
                         chr_lens = med_chr_lens)
    
    return(out)
  })
  names(bin_list) = names(gwas_true[[COUNTER_1]])
  
  return(bin_list)
})
names(plot_dat) = names(gwas_true)
```

Read in significance levels from permutations


```r
perm_file = here::here("data/20211109_permutation_mins.csv")
perms_df_mins = readr::read_csv(perm_file,
                                col_types = c("ccid"))
```

Plot

```r
plot_dir = here::here(file.path("book/plots/20211109/gwls_results", filter_type))
# Plot
lapply(seq_along(plot_dat), function(COUNTER_PHENO){
  lapply(seq_along(plot_dat[[COUNTER_PHENO]]), function(COUNTER_BINL){
    # Get pheno
    PHENO = names(plot_dat)[COUNTER_PHENO]
    # Get bin length
    BINL = names(plot_dat[[COUNTER_PHENO]])[COUNTER_BINL] %>% 
      as.numeric()
    # Get significant level
    SIG_LEVEL = perms_df_mins %>%
      dplyr::filter(SITE_FILTER == filter_type, TARGET_PHENO == PHENO, BIN_LENGTH == BINL) %>% 
      dplyr::pull(MIN_P) %>% 
      -log10(.)
    
    # Plot
    out = plot_man(plot_dat[[COUNTER_PHENO]][[COUNTER_BINL]],
             phenotype = PHENO,
             bin_length = BINL,
             gwas_pal = gwas_pal[2:3],
             med_chr_lens = med_chr_lens,
             sig_line = SIG_LEVEL) +
      ylim(0,7)
    
    ggsave(file.path(plot_dir, paste(PHENO, "_", BINL, ".png", sep = "")),
           out,
           device = "png",
           width = 9.6,
           height = 6,
           units = "in",
           dpi = 400)
    
    out
  })
})
#> [[1]]
#> [[1]][[1]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-1.png" width="672" />

```
#> 
#> [[1]][[2]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-2.png" width="672" />

```
#> 
#> [[1]][[3]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-3.png" width="672" />

```
#> 
#> [[1]][[4]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-4.png" width="672" />

```
#> 
#> 
#> [[2]]
#> [[2]][[1]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-5.png" width="672" />

```
#> 
#> [[2]][[2]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-6.png" width="672" />

```
#> 
#> [[2]][[3]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-7.png" width="672" />

```
#> 
#> [[2]][[4]]
```

<img src="05-GWLS_results_files/figure-html/unnamed-chunk-9-8.png" width="672" />
