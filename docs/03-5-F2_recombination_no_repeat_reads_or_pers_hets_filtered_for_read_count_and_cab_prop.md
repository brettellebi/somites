---
zotero: "PhD"
---

# F2 recombination blocks

Exclusions:
-   reads overlapping *HdrR* repeat regions
-   regions of persistent heterozygosity in the MIKK panel
-   filtered based on read count and proportion of Cab)

Snakefile for generating figures: https://github.com/brettellebi/somites/blob/master/workflow/rules/05_F2_recomb_blocks.smk


```r
library(here)
#> here() starts at /Users/brettell/Documents/Repositories/somites
site_filter = "no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop"
```

## Base coverage

### Total


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/base_cov_total.png"))
```

<img src="/Users/brettell/Documents/Repositories/somites/book/plots/snakemake/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/5000/base_cov_total.png" width="100%" />

### By chromosome


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/base_cov_by_chrom.png"))
```

<img src="/Users/brettell/Documents/Repositories/somites/book/plots/snakemake/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/5000/base_cov_by_chrom.png" width="100%" />

## Proportion of sites

### Total


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/prop_sites_total.png"))
```

<img src="/Users/brettell/Documents/Repositories/somites/book/plots/snakemake/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/5000/prop_sites_total.png" width="100%" />

### By chromosome


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/prop_sites_by_chrom.png"))
```

<img src="/Users/brettell/Documents/Repositories/somites/book/plots/snakemake/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/5000/prop_sites_by_chrom.png" width="100%" />

## Karyoplots

### No missing blocks


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/karyoplot_no_missing.png"))
```

<img src="/Users/brettell/Documents/Repositories/somites/book/plots/snakemake/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/5000/karyoplot_no_missing.png" width="100%" />

### With missing blocks


```r
knitr::include_graphics(here::here("book/plots/snakemake", site_filter, "5000/karyoplot_with_missing.png"))
```

<img src="/Users/brettell/Documents/Repositories/somites/book/plots/snakemake/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/5000/karyoplot_with_missing.png" width="100%" />


