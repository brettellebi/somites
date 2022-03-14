---
zotero: PhD
---

# GWAS results


```r
library(here)
source(here::here("book/source/04-Association_testing.R"))
```

## Notes

* 20211104 association test performed on full sites file. Results here: `/nfs/research/birney/users/ian/somites/association_testing/20211104_true/results`
* 20211109 association test performed on sites excluding those overlapping repeat regions. Results here: `/nfs/research/birney/users/ian/somites/association_testing/20211109_true/results`
* 20220214 association tests performed on all filter types and including PSM size as a third phenotype.

## Snakemake rules

Snakemake rules for running the GWAS over these phenotypes: <https://github.com/brettellebi/somites/blob/master/workflow/rules/07_assocation_testing.smk>

## Results


```r
DATE_OF_ASSOC_TEST = 20220214
```

### All sites


```r
SITE_FILTER = "all_sites"
INVERSE_NORM = "TRUE"
```

#### Intercept


```r
TARGET_PHENO = "intercept"
```


```r
COVARIATES = "None"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/intercept/None/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/intercept/None/TRUE/20000.png" width="100%" />


```r
COVARIATES = "Microscope"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/intercept/Microscope/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/intercept/Microscope/TRUE/20000.png" width="100%" />

#### Mean


```r
TARGET_PHENO = "mean"
```


```r
COVARIATES = "None"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/mean/None/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/mean/None/TRUE/20000.png" width="100%" />


```r
COVARIATES = "Microscope"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/mean/Microscope/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/mean/Microscope/TRUE/20000.png" width="100%" />

#### PSM


```r
TARGET_PHENO = "unsegmented_psm_area"
```


```r
COVARIATES = "None"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/unsegmented_psm_area/None/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/unsegmented_psm_area/None/TRUE/20000.png" width="100%" />


```r
COVARIATES = "Microscope"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/unsegmented_psm_area/Microscope/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/all_sites/unsegmented_psm_area/Microscope/TRUE/20000.png" width="100%" />

### All filters


```r
SITE_FILTER = "no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop"
INVERSE_NORM = "TRUE"
```

#### Intercept


```r
TARGET_PHENO = "intercept"
```


```r
COVARIATES = "None"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/intercept/None/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/intercept/None/TRUE/20000.png" width="100%" />


```r
COVARIATES = "Microscope"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/intercept/Microscope/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/intercept/Microscope/TRUE/20000.png" width="100%" />


#### Mean


```r
TARGET_PHENO = "mean"
```


```r
COVARIATES = "None"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/mean/None/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/mean/None/TRUE/20000.png" width="100%" />


```r
COVARIATES = "Microscope"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/mean/Microscope/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/mean/Microscope/TRUE/20000.png" width="100%" />


#### PSM


```r
TARGET_PHENO = "unsegmented_psm_area"
```


```r
COVARIATES = "None"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/unsegmented_psm_area/None/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/unsegmented_psm_area/None/TRUE/20000.png" width="100%" />


```r
COVARIATES = "Microscope"
```


```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "5000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/unsegmented_psm_area/Microscope/TRUE/5000.png" width="100%" />

```r
knitr::include_graphics(here::here("book/plots/manhattans", DATE_OF_ASSOC_TEST, SITE_FILTER, TARGET_PHENO, COVARIATES, INVERSE_NORM, "20000.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/manhattans/20220214/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/unsegmented_psm_area/Microscope/TRUE/20000.png" width="100%" />


