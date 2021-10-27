---
zotero: PhD
---

# Association testing


```r
library(here)
source(here::here("book/source/04-Association_testing.R"))
```

## Phenotype data

First 400 phenotype data: <https://github.com/brettellebi/somites/tree/master/data/20210917_First400_F2_DF.xlsx>.

>I also attach here the DataFrame of the phenotyping of the first 400 F2s (First400_F2_DF.xlsx)...

>Of interest for us for the association testing are the:
 intercept period -> intercept in the table
 mean period -> mean in the table


```r
pheno_file = here::here("data/20210917_First400_F2_DF.xlsx")
pheno = readxl::read_xlsx(pheno_file)

# Which samples are missing phenotypes? 
which(!1:400 %in% as.integer(pheno$fish %>% str_remove("KC")))
#> [1]  23 278
```

