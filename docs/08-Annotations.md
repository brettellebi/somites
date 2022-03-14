# Annotations of GWAS hits


```r
library(here)
#> here() starts at /hps/software/users/birney/ian/repos/somites
library(tidyverse)
#> ── Attaching packages ─────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.3     ✓ dplyr   1.0.7
#> ✓ tidyr   1.1.3     ✓ stringr 1.4.0
#> ✓ readr   2.0.0     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
library(biomaRt)
```


```r
DATE_OF_ASSOC_TEST = "20220214"
SITE_FILTER = "all_sites"
COVARIATES = "None"
INVERSE_NORM = "TRUE"
BIN_LENGTH = 5000
OUT_DIR = here::here("data", DATE_OF_ASSOC_TEST, "annotations")
```

## Intercept


```r
TARGET_PHENO = "intercept"
```

### Read in data


```r
GWAS_RESULTS = file.path("/hps/nobackup/birney/users/ian/somites/association_testing",
                         DATE_OF_ASSOC_TEST,
                         SITE_FILTER,
                         "true_results",
                         TARGET_PHENO,
                         COVARIATES,
                         INVERSE_NORM,
                         paste(BIN_LENGTH, ".rds", sep = ""))

RESULTS = readRDS(GWAS_RESULTS)

PERMS_PATHS = file.path("/hps/nobackup/birney/users/ian/somites/association_testing/",
                          DATE_OF_ASSOC_TEST,
                          SITE_FILTER,
                          "permutations",
                          TARGET_PHENO,
                          COVARIATES,
                          INVERSE_NORM,
                          BIN_LENGTH)
PERMS = list.files(PERMS_PATHS, full.names = T)
names(PERMS) = PERMS %>% 
  basename %>% 
  stringr::str_remove(".rds")

if ("p_value_REML" %in% colnames(RESULTS$results)){
  P_COL = "p_value_REML"
} else {
  P_COL = "p_value_REML.1"
}

# Rename column in results

if (P_COL == "p_value_REML.1"){
  RESULTS$results = RESULTS$results %>% 
    dplyr::rename(p_value_REML = p_value_REML.1)
}

# Read in permutation results and get `SIG_LEVEL`

PERM_LIST = purrr::map(PERMS, function(PERM){
  readRDS(PERM)
})

perm_df = purrr::map_dfr(PERM_LIST, function(PERM){
  OUT = tibble::tibble(MIN_P = PERM$results %>% 
                         dplyr::select(dplyr::all_of(P_COL)) %>%
                         min(., na.rm = T)
  )
}, .id = "SEED")

# Get minimum
SIG_LEVEL = min(perm_df$MIN_P)
```

### Pull significant loci


```r
SIG_LOCS = RESULTS$results %>% 
  dplyr::filter(p_value_REML < SIG_LEVEL) %>% 
  dplyr::select(CHROM = Chr,
                BIN_START = pos) %>% 
  dplyr::mutate(BIN_END = BIN_START + BIN_LENGTH - 1)
                
```


### Get annotations


```r
## Select dataset
olat_mart = biomaRt::useEnsembl(biomart = "ensembl", 
                                dataset = "olatipes_gene_ensembl", 
                                mirror = "uswest")

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

# convert hits to genomic ranges
sig_loc_r = SIG_LOCS %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM",
                                          start.field = "BIN_START",
                                          end.field = "BIN_END",
                                          ignore.strand = T)


# find overlaps
olaps = GenomicRanges::findOverlaps(sig_loc_r, olat_genes_r)

# Pull out data frame of hits
hits = olat_genes[unique(olaps@to), ]

hits %>% 
  DT::datatable(.)
```

```{=html}
<div id="htmlwidget-37e562e844253ac993b1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-37e562e844253ac993b1">{"x":{"filter":"none","data":[["52786"],["3"],[34635715],[34636717],["ENSORLG00000026609"],[""],["ENSORLE00000271723"],[""],[1],[34635715],[34636717]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>chromosome_name<\/th>\n      <th>start_position<\/th>\n      <th>end_position<\/th>\n      <th>ensembl_gene_id<\/th>\n      <th>hgnc_symbol<\/th>\n      <th>ensembl_exon_id<\/th>\n      <th>description<\/th>\n      <th>strand<\/th>\n      <th>transcript_start<\/th>\n      <th>transcript_end<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
# Save to file
readr::write_csv(hits,
                 file.path(OUT_DIR, paste(TARGET_PHENO, ".csv", sep = "")))
```

## PSM area


```r
TARGET_PHENO = "unsegmented_psm_area"
```

### Read in data


```r
GWAS_RESULTS = file.path("/hps/nobackup/birney/users/ian/somites/association_testing",
                         DATE_OF_ASSOC_TEST,
                         SITE_FILTER,
                         "true_results",
                         TARGET_PHENO,
                         COVARIATES,
                         INVERSE_NORM,
                         paste(BIN_LENGTH, ".rds", sep = ""))

RESULTS = readRDS(GWAS_RESULTS)

PERMS_PATHS = file.path("/hps/nobackup/birney/users/ian/somites/association_testing/",
                          DATE_OF_ASSOC_TEST,
                          SITE_FILTER,
                          "permutations",
                          TARGET_PHENO,
                          COVARIATES,
                          INVERSE_NORM,
                          BIN_LENGTH)
PERMS = list.files(PERMS_PATHS, full.names = T)
names(PERMS) = PERMS %>% 
  basename %>% 
  stringr::str_remove(".rds")

if ("p_value_REML" %in% colnames(RESULTS$results)){
  P_COL = "p_value_REML"
} else {
  P_COL = "p_value_REML.1"
}

# Rename column in results

if (P_COL == "p_value_REML.1"){
  RESULTS$results = RESULTS$results %>% 
    dplyr::rename(p_value_REML = p_value_REML.1)
}

# Read in permutation results and get `SIG_LEVEL`

PERM_LIST = purrr::map(PERMS, function(PERM){
  readRDS(PERM)
})

perm_df = purrr::map_dfr(PERM_LIST, function(PERM){
  OUT = tibble::tibble(MIN_P = PERM$results %>% 
                         dplyr::select(dplyr::all_of(P_COL)) %>%
                         min(., na.rm = T)
  )
}, .id = "SEED")

# Get minimum
SIG_LEVEL = min(perm_df$MIN_P)
```

### Pull significant loci


```r
SIG_LOCS = RESULTS$results %>% 
  dplyr::filter(p_value_REML < SIG_LEVEL) %>% 
  dplyr::select(CHROM = Chr,
                BIN_START = pos) %>% 
  dplyr::mutate(BIN_END = BIN_START + BIN_LENGTH - 1)
                
```


### Get annotations


```r
# convert hits to genomic ranges
sig_loc_r = SIG_LOCS %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM",
                                          start.field = "BIN_START",
                                          end.field = "BIN_END",
                                          ignore.strand = T)


# find overlaps
olaps = GenomicRanges::findOverlaps(sig_loc_r, olat_genes_r)

# Pull out data frame of hits
hits = olat_genes[unique(olaps@to), ]

hits %>% 
  DT::datatable(.)
```

```{=html}
<div id="htmlwidget-a965ddca991e0e9fe965" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-a965ddca991e0e9fe965">{"x":{"filter":"none","data":[["106182","106183","106184","106185","106186","106187","106188","106189","106190","106191","106192","106193","106194","106195","106196"],["3","3","3","3","3","3","3","3","3","3","3","3","3","3","3"],[22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087,22963087],[22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992674],["ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678","ENSORLG00000009678"],["","","","","","","","","","","","","","",""],["ENSORLE00000270694","ENSORLE00000109779","ENSORLE00000109778","ENSORLE00000109768","ENSORLE00000109762","ENSORLE00000269835","ENSORLE00000302288","ENSORLE00000249511","ENSORLE00000109760","ENSORLE00000109757","ENSORLE00000109761","ENSORLE00000138828","ENSORLE00000109778","ENSORLE00000109768","ENSORLE00000109762"],["synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]","synaptotagmin-9 [Source:NCBI gene;Acc:101174010]"],[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],[22963087,22963087,22963087,22963087,22963087,22963087,22963087,22968889,22968889,22968889,22968889,22968889,22968889,22968889,22968889],[22992674,22992674,22992674,22992674,22992674,22992674,22992674,22992573,22992573,22992573,22992573,22992573,22992573,22992573,22992573]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>chromosome_name<\/th>\n      <th>start_position<\/th>\n      <th>end_position<\/th>\n      <th>ensembl_gene_id<\/th>\n      <th>hgnc_symbol<\/th>\n      <th>ensembl_exon_id<\/th>\n      <th>description<\/th>\n      <th>strand<\/th>\n      <th>transcript_start<\/th>\n      <th>transcript_end<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
# Save to file
readr::write_csv(hits,
                 file.path(OUT_DIR, paste(TARGET_PHENO, ".csv", sep = "")))
```

