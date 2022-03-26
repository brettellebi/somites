# Annotations of GWAS hits


```r
library(here)
library(tidyverse)
library(biomaRt)
```


```r
DATE_OF_ASSOC_TEST = "20220321"
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
<div id="htmlwidget-bc517b9517d1f250c76e" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bc517b9517d1f250c76e">{"x":{"filter":"none","vertical":false,"data":[["56340","56341","56342","56343","56344","56345","56346","56347","56348"],["19","19","19","19","19","19","19","19","19"],[25403379,25403379,25403379,25403379,25403379,25403379,25403379,25403379,25403379],[25426536,25426536,25426536,25426536,25426536,25426536,25426536,25426536,25426536],["ENSORLG00000030311","ENSORLG00000030311","ENSORLG00000030311","ENSORLG00000030311","ENSORLG00000030311","ENSORLG00000030311","ENSORLG00000030311","ENSORLG00000030311","ENSORLG00000030311"],["","","","","","","","",""],["ENSORLE00000214856","ENSORLE00000214853","ENSORLE00000214851","ENSORLE00000261730","ENSORLE00000237352","ENSORLE00000272226","ENSORLE00000228790","ENSORLE00000271898","ENSORLE00000308291"],["","","","","","","","",""],[1,1,1,1,1,1,1,1,1],[25403379,25403379,25403379,25403379,25403379,25403379,25403379,25403379,25403379],[25426536,25426536,25426536,25426536,25426536,25426536,25426536,25426536,25426536]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>chromosome_name<\/th>\n      <th>start_position<\/th>\n      <th>end_position<\/th>\n      <th>ensembl_gene_id<\/th>\n      <th>hgnc_symbol<\/th>\n      <th>ensembl_exon_id<\/th>\n      <th>description<\/th>\n      <th>strand<\/th>\n      <th>transcript_start<\/th>\n      <th>transcript_end<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
# Save to file
readr::write_csv(hits,
                 file.path(OUT_DIR, paste(TARGET_PHENO, ".csv", sep = "")))
```


