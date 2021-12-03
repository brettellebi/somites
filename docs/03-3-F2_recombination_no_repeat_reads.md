---
zotero: "PhD"
---

# F2 recombination blocks (excluding READS overlapping HdrR repeat regions)

Snakefile for generating figures: https://github.com/brettellebi/somites/blob/master/workflow/rules/05_F2_recomb_blocks.smk


```r
library(here)
#> here() starts at /hps/software/users/birney/ian/repos/somites
site_filter = "no_repeat_reads"
```

## Base coverage

### Total

#### Bin size: 5000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/base_cov_total.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/5000/base_cov_total.png" width="2000" />

#### Bin size: 20000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "20000/base_cov_total.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/20000/base_cov_total.png" width="2000" />

### By chromosome

#### Bin size: 5000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/base_cov_by_chrom.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/5000/base_cov_by_chrom.png" width="3200" />

#### Bin size: 20000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "20000/base_cov_by_chrom.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/20000/base_cov_by_chrom.png" width="3200" />

## Proportion of sites

### Total

#### Bin size: 5000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/prop_sites_total.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/5000/prop_sites_total.png" width="2000" />

#### Bin size: 20000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "20000/prop_sites_total.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/20000/prop_sites_total.png" width="2000" />

### By chromosome

#### Bin size: 5000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/prop_sites_by_chrom.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/5000/prop_sites_by_chrom.png" width="3200" />

#### Bin size: 20000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "20000/prop_sites_by_chrom.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/20000/prop_sites_by_chrom.png" width="3200" />

## Karyoplots

### No missing blocks

#### Bin size: 5000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/karyoplot_no_missing.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/5000/karyoplot_no_missing.png" width="3900" />

#### Bin size: 20000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "20000/karyoplot_no_missing.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/20000/karyoplot_no_missing.png" width="3900" />

### With missing blocks

#### Bin size: 5000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/karyoplot_with_missing.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/5000/karyoplot_with_missing.png" width="3900" />

#### Bin size: 20000 bp


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "20000/karyoplot_with_missing.png"))
```

<img src="/hps/software/users/birney/ian/repos/somites/book/plots/snakemake/no_repeat_reads/20000/karyoplot_with_missing.png" width="3900" />

