---
zotero: PhD
---

# Association testing


```r
library(here)
source(here::here("book/source/04-Association_testing.R"))
```

## Read in F2 recombination blocks

### Read data


```r
in_dir = "/nfs/research/birney/users/ian/somites/recombination_blocks/20211027"

in_files = list.files(in_dir, pattern = "F2_", full.names = T)

# Read into list
data_list = purrr::map(in_files, function(FILE){
  out = readr::read_tsv(FILE,
                        col_types = "ciiidii")
})
# Set names as bin length
names(data_list) = basename(in_files) %>% 
  stringr::str_split("_", simplify = T) %>% 
  subset(select = 2) %>% 
  stringr::str_remove(".txt")
# Reorder
data_list = data_list[order(as.numeric(names(data_list)))]

counter = 0
df_list = purrr::map(data_list, function(data){
  counter <<- counter + 1
  # set bin length
  bin_length = as.numeric(names(data_list)[counter])
  # add bin start and end coordinates
  df = data %>% 
    dplyr::mutate(SAMPLE = basename(sample) %>% 
                    stringr::str_remove(".txt") %>% 
                    as.numeric(.),
                  BIN_START = (bin - 1) * bin_length + 1,
                  BIN_END = bin * bin_length) %>% 
    # recode state to make 0 == "Cab"
    dplyr::mutate(STATE = dplyr::recode(state,
                                        `0` = 2,
                                        `1` = 1,
                                        `2` = 0)) %>% 
    dplyr::select(SAMPLE, CHROM = chr, BIN = bin, BIN_START, BIN_END, STATE)
  
  return(df)
})

```

### How many possible blocks have at least one call? 

#### Count total existing bins


```r
bin_lengths = as.integer(names(df_list))
names(bin_lengths) = bin_lengths

n_bins = purrr::map(bin_lengths, function(BIN_LENGTH){
  purrr::map_int(med_chr_lens$end, function(CHR_END){
    out = seq(from = 1, to = CHR_END, by = BIN_LENGTH)
    
    return(length(out))
  })
})

n_bins_total = purrr::map_int(n_bins, sum)
```

#### Proportion of total bins with calls


```r
# Get total number of samples
n_samples = df_list$`5000`$SAMPLE %>%
  unique(.) %>% 
  length(.)

# Build DF of bins
n_bins_df = purrr::map_dfr(1:length(df_list), function(COUNTER){
  # Bin length
  bin_length = as.numeric(names(df_list)[COUNTER])
  # Number of total bins
  n_bins = n_bins_total[COUNTER]
  # Number of bins with calls
  n_bins_with_calls = df_list[[COUNTER]] %>% 
    dplyr::distinct(CHROM, BIN_START) %>% 
    nrow(.)
  # Number of bins with calls for each 
  n_bins_no_missing = df_list[[COUNTER]] %>% 
    dplyr::count(CHROM, BIN_START) %>%
    dplyr::filter(n == n_samples) %>%
    nrow(.)
  
  # Build final data frame
  out = tibble::tibble(BIN_LENGTH = bin_length,
                       N_BINS = n_bins,
                       N_BINS_WITH_CALLS = n_bins_with_calls,
                       N_BINS_NO_MISSING = n_bins_no_missing) %>% 
    dplyr::mutate(PROP_BINS_WITH_CALLS = N_BINS_WITH_CALLS / N_BINS,
                  PROP_BINS_NO_MISSING = N_BINS_NO_MISSING / N_BINS )
  
  return(out)
})

DT::datatable(n_bins_df)
```

```{=html}
<div id="htmlwidget-caa78431e141ff308e71" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-caa78431e141ff308e71">{"x":{"filter":"none","data":[["1","2","3","4"],[5000,10000,15000,20000],[146819,73414,48947,36712],[111237,62489,43734,33761],[10256,11886,11909,11520],[0.757647171006477,0.851186422208298,0.893497047827242,0.919617563739377],[0.0698547190758689,0.161903724085324,0.243303981857928,0.313793854870342]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>BIN_LENGTH<\/th>\n      <th>N_BINS<\/th>\n      <th>N_BINS_WITH_CALLS<\/th>\n      <th>N_BINS_NO_MISSING<\/th>\n      <th>PROP_BINS_WITH_CALLS<\/th>\n      <th>PROP_BINS_NO_MISSING<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


### Merge recombination blocks

#### List with final recoded genotypes (including NAs)


```r
gt_list = purrr::map(df_list, function(BIN_LENGTH_DF){
  
  # Widen data frame
  gt_final = BIN_LENGTH_DF %>% 
    tidyr::pivot_wider(names_from = SAMPLE, values_from = STATE)
  
  # Pull out matrix of genotypes
  gt_mat = as.matrix(gt_final[, 5:ncol(gt_final)])
  
  # Get indexes of loci with > 1 genotype
  bins_to_keep = logical()
  
  for (ROW in 1:nrow(gt_mat)){
    # get unique values in each row
    out = unique(gt_mat[ROW, ])
    # remove NAs
    out = out[!is.na(out)]
    # if more than one value, return TRUE
    if (length(out) > 1) {
      bins_to_keep[ROW] = TRUE
    }
    # if just one value (i.e. if all samples are the same genotype at that locus), return false 
    else {
      bins_to_keep[ROW] = FALSE
    }
  }
  
  # filter gt_final
  gt_filt = gt_final %>% 
    dplyr::filter(bins_to_keep) %>% 
    # recode genotypes to -1, 0, 1
    dplyr::mutate(dplyr::across(-c("CHROM", "BIN", "BIN_START", "BIN_END"),
                                ~dplyr::recode(.x,
                                               `0` = -1,
                                               `1` = 0,
                                               `2` = 1))) %>% 
    # order
    dplyr::arrange(CHROM, BIN_START)  
  
  return(gt_filt)
}) 


gt_list$`20000` %>% 
  head(.) %>% 
  DT::datatable(.) 
```

```{=html}
<div id="htmlwidget-b20a1e645f0095da9596" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-b20a1e645f0095da9596">{"x":{"filter":"none","data":[["1","2","3","4","5","6"],[1,1,1,1,1,1],[1,2,3,4,5,6],[1,20001,40001,60001,80001,100001],[20000,40000,60000,80000,100000,120000],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,null,0,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[-1,-1,-1,-1,-1,-1],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,null,null,-1,-1,-1],[0,-1,-1,-1,-1,-1],[null,1,null,1,1,1],[-1,null,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[1,1,null,1,1,1],[-1,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[-1,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,0,-1,-1,-1],[0,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[null,null,0,1,1,1],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[1,1,1,1,1,1],[1,1,null,1,1,1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,-1,-1,-1,-1],[0,-1,-1,-1,-1,-1],[0,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[-1,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[1,1,1,1,1,1],[-1,-1,null,-1,-1,-1],[-1,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[-1,-1,null,-1,-1,-1],[0,-1,-1,-1,-1,-1],[0,0,null,0,0,0],[0,0,0,0,0,0],[0,null,null,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[1,1,1,1,1,1],[0,0,null,0,0,0],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,0,0,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,-1,-1,-1,-1],[0,0,0,0,0,0],[1,1,null,1,1,1],[-1,-1,-1,-1,-1,-1],[0,0,0,-1,-1,-1],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,0,0,0,0,0],[null,0,0,1,1,1],[0,-1,-1,-1,-1,-1],[1,1,null,1,1,1],[0,0,null,0,0,0],[0,-1,-1,-1,-1,-1],[0,0,0,0,0,0],[0,0,0,0,0,0],[-1,-1,-1,-1,-1,-1],[0,0,0,0,0,0],[0,-1,-1,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,0,1,1,1],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,-1,-1,-1,-1],[0,0,0,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,null,0,0,0,0],[-1,-1,-1,-1,-1,-1],[0,0,null,0,0,0],[0,0,0,0,0,0],[0,0,null,0,0,0],[null,0,0,0,0,0],[null,0,0,1,null,1],[0,-1,-1,-1,-1,-1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,0,0,0,0],[0,0,null,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[-1,-1,-1,-1,-1,-1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,null,0,0,0],[null,null,0,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[1,1,null,1,1,1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[null,null,null,null,null,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,0,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[1,1,null,1,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[1,1,null,1,1,1],[-1,-1,null,-1,-1,-1],[1,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,null,-1],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[1,1,null,1,1,1],[1,1,null,1,1,1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[1,1,null,1,1,1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[1,1,null,1,1,1],[-1,-1,null,-1,-1,-1],[1,1,null,1,1,1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[1,1,null,1,1,1],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,null],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[1,1,null,1,1,1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[1,1,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[1,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,1,null,1,1,1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[1,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,null],[0,0,null,0,0,null],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[1,1,null,1,1,1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,-1,null,-1,-1,-1],[0,0,null,0,0,0],[0,0,null,0,0,0],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1],[0,0,null,0,0,0],[-1,-1,null,-1,-1,-1]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>CHROM<\/th>\n      <th>BIN<\/th>\n      <th>BIN_START<\/th>\n      <th>BIN_END<\/th>\n      <th>1<\/th>\n      <th>2<\/th>\n      <th>4<\/th>\n      <th>6<\/th>\n      <th>8<\/th>\n      <th>9<\/th>\n      <th>10<\/th>\n      <th>11<\/th>\n      <th>12<\/th>\n      <th>13<\/th>\n      <th>14<\/th>\n      <th>15<\/th>\n      <th>16<\/th>\n      <th>17<\/th>\n      <th>18<\/th>\n      <th>20<\/th>\n      <th>21<\/th>\n      <th>22<\/th>\n      <th>23<\/th>\n      <th>24<\/th>\n      <th>25<\/th>\n      <th>26<\/th>\n      <th>27<\/th>\n      <th>29<\/th>\n      <th>30<\/th>\n      <th>31<\/th>\n      <th>32<\/th>\n      <th>33<\/th>\n      <th>34<\/th>\n      <th>35<\/th>\n      <th>36<\/th>\n      <th>37<\/th>\n      <th>38<\/th>\n      <th>40<\/th>\n      <th>41<\/th>\n      <th>42<\/th>\n      <th>43<\/th>\n      <th>44<\/th>\n      <th>45<\/th>\n      <th>46<\/th>\n      <th>47<\/th>\n      <th>48<\/th>\n      <th>49<\/th>\n      <th>50<\/th>\n      <th>51<\/th>\n      <th>53<\/th>\n      <th>54<\/th>\n      <th>55<\/th>\n      <th>57<\/th>\n      <th>58<\/th>\n      <th>59<\/th>\n      <th>60<\/th>\n      <th>61<\/th>\n      <th>62<\/th>\n      <th>63<\/th>\n      <th>64<\/th>\n      <th>65<\/th>\n      <th>66<\/th>\n      <th>67<\/th>\n      <th>68<\/th>\n      <th>69<\/th>\n      <th>70<\/th>\n      <th>71<\/th>\n      <th>72<\/th>\n      <th>73<\/th>\n      <th>74<\/th>\n      <th>75<\/th>\n      <th>76<\/th>\n      <th>77<\/th>\n      <th>78<\/th>\n      <th>79<\/th>\n      <th>80<\/th>\n      <th>81<\/th>\n      <th>82<\/th>\n      <th>83<\/th>\n      <th>84<\/th>\n      <th>85<\/th>\n      <th>86<\/th>\n      <th>87<\/th>\n      <th>88<\/th>\n      <th>89<\/th>\n      <th>90<\/th>\n      <th>91<\/th>\n      <th>92<\/th>\n      <th>93<\/th>\n      <th>94<\/th>\n      <th>95<\/th>\n      <th>96<\/th>\n      <th>97<\/th>\n      <th>98<\/th>\n      <th>99<\/th>\n      <th>100<\/th>\n      <th>101<\/th>\n      <th>102<\/th>\n      <th>103<\/th>\n      <th>104<\/th>\n      <th>105<\/th>\n      <th>106<\/th>\n      <th>107<\/th>\n      <th>108<\/th>\n      <th>109<\/th>\n      <th>110<\/th>\n      <th>111<\/th>\n      <th>112<\/th>\n      <th>113<\/th>\n      <th>114<\/th>\n      <th>115<\/th>\n      <th>116<\/th>\n      <th>117<\/th>\n      <th>118<\/th>\n      <th>119<\/th>\n      <th>120<\/th>\n      <th>122<\/th>\n      <th>123<\/th>\n      <th>124<\/th>\n      <th>125<\/th>\n      <th>126<\/th>\n      <th>127<\/th>\n      <th>128<\/th>\n      <th>129<\/th>\n      <th>130<\/th>\n      <th>131<\/th>\n      <th>132<\/th>\n      <th>133<\/th>\n      <th>135<\/th>\n      <th>136<\/th>\n      <th>137<\/th>\n      <th>138<\/th>\n      <th>139<\/th>\n      <th>140<\/th>\n      <th>141<\/th>\n      <th>142<\/th>\n      <th>143<\/th>\n      <th>144<\/th>\n      <th>145<\/th>\n      <th>146<\/th>\n      <th>147<\/th>\n      <th>148<\/th>\n      <th>149<\/th>\n      <th>150<\/th>\n      <th>151<\/th>\n      <th>152<\/th>\n      <th>153<\/th>\n      <th>154<\/th>\n      <th>156<\/th>\n      <th>157<\/th>\n      <th>159<\/th>\n      <th>161<\/th>\n      <th>162<\/th>\n      <th>163<\/th>\n      <th>164<\/th>\n      <th>165<\/th>\n      <th>166<\/th>\n      <th>167<\/th>\n      <th>168<\/th>\n      <th>169<\/th>\n      <th>170<\/th>\n      <th>171<\/th>\n      <th>172<\/th>\n      <th>174<\/th>\n      <th>175<\/th>\n      <th>176<\/th>\n      <th>178<\/th>\n      <th>179<\/th>\n      <th>180<\/th>\n      <th>181<\/th>\n      <th>182<\/th>\n      <th>183<\/th>\n      <th>184<\/th>\n      <th>185<\/th>\n      <th>186<\/th>\n      <th>187<\/th>\n      <th>188<\/th>\n      <th>189<\/th>\n      <th>190<\/th>\n      <th>191<\/th>\n      <th>192<\/th>\n      <th>193<\/th>\n      <th>194<\/th>\n      <th>195<\/th>\n      <th>196<\/th>\n      <th>197<\/th>\n      <th>198<\/th>\n      <th>199<\/th>\n      <th>200<\/th>\n      <th>201<\/th>\n      <th>202<\/th>\n      <th>203<\/th>\n      <th>204<\/th>\n      <th>205<\/th>\n      <th>206<\/th>\n      <th>207<\/th>\n      <th>208<\/th>\n      <th>209<\/th>\n      <th>210<\/th>\n      <th>212<\/th>\n      <th>213<\/th>\n      <th>214<\/th>\n      <th>215<\/th>\n      <th>216<\/th>\n      <th>218<\/th>\n      <th>219<\/th>\n      <th>220<\/th>\n      <th>221<\/th>\n      <th>222<\/th>\n      <th>223<\/th>\n      <th>224<\/th>\n      <th>225<\/th>\n      <th>226<\/th>\n      <th>227<\/th>\n      <th>228<\/th>\n      <th>229<\/th>\n      <th>230<\/th>\n      <th>231<\/th>\n      <th>232<\/th>\n      <th>233<\/th>\n      <th>234<\/th>\n      <th>235<\/th>\n      <th>236<\/th>\n      <th>237<\/th>\n      <th>238<\/th>\n      <th>239<\/th>\n      <th>240<\/th>\n      <th>241<\/th>\n      <th>242<\/th>\n      <th>243<\/th>\n      <th>244<\/th>\n      <th>245<\/th>\n      <th>247<\/th>\n      <th>248<\/th>\n      <th>249<\/th>\n      <th>250<\/th>\n      <th>251<\/th>\n      <th>252<\/th>\n      <th>253<\/th>\n      <th>254<\/th>\n      <th>255<\/th>\n      <th>256<\/th>\n      <th>259<\/th>\n      <th>260<\/th>\n      <th>261<\/th>\n      <th>262<\/th>\n      <th>263<\/th>\n      <th>264<\/th>\n      <th>265<\/th>\n      <th>266<\/th>\n      <th>267<\/th>\n      <th>268<\/th>\n      <th>269<\/th>\n      <th>270<\/th>\n      <th>271<\/th>\n      <th>272<\/th>\n      <th>273<\/th>\n      <th>274<\/th>\n      <th>275<\/th>\n      <th>276<\/th>\n      <th>277<\/th>\n      <th>278<\/th>\n      <th>280<\/th>\n      <th>281<\/th>\n      <th>282<\/th>\n      <th>283<\/th>\n      <th>284<\/th>\n      <th>285<\/th>\n      <th>286<\/th>\n      <th>287<\/th>\n      <th>288<\/th>\n      <th>289<\/th>\n      <th>290<\/th>\n      <th>291<\/th>\n      <th>292<\/th>\n      <th>293<\/th>\n      <th>294<\/th>\n      <th>295<\/th>\n      <th>296<\/th>\n      <th>297<\/th>\n      <th>298<\/th>\n      <th>299<\/th>\n      <th>300<\/th>\n      <th>301<\/th>\n      <th>302<\/th>\n      <th>303<\/th>\n      <th>304<\/th>\n      <th>306<\/th>\n      <th>307<\/th>\n      <th>308<\/th>\n      <th>309<\/th>\n      <th>310<\/th>\n      <th>311<\/th>\n      <th>312<\/th>\n      <th>313<\/th>\n      <th>314<\/th>\n      <th>315<\/th>\n      <th>316<\/th>\n      <th>317<\/th>\n      <th>318<\/th>\n      <th>319<\/th>\n      <th>320<\/th>\n      <th>321<\/th>\n      <th>322<\/th>\n      <th>323<\/th>\n      <th>324<\/th>\n      <th>326<\/th>\n      <th>327<\/th>\n      <th>328<\/th>\n      <th>329<\/th>\n      <th>330<\/th>\n      <th>331<\/th>\n      <th>332<\/th>\n      <th>333<\/th>\n      <th>334<\/th>\n      <th>335<\/th>\n      <th>336<\/th>\n      <th>337<\/th>\n      <th>338<\/th>\n      <th>339<\/th>\n      <th>340<\/th>\n      <th>341<\/th>\n      <th>342<\/th>\n      <th>343<\/th>\n      <th>344<\/th>\n      <th>345<\/th>\n      <th>346<\/th>\n      <th>347<\/th>\n      <th>348<\/th>\n      <th>349<\/th>\n      <th>350<\/th>\n      <th>351<\/th>\n      <th>352<\/th>\n      <th>353<\/th>\n      <th>354<\/th>\n      <th>355<\/th>\n      <th>356<\/th>\n      <th>357<\/th>\n      <th>358<\/th>\n      <th>359<\/th>\n      <th>360<\/th>\n      <th>361<\/th>\n      <th>362<\/th>\n      <th>363<\/th>\n      <th>364<\/th>\n      <th>365<\/th>\n      <th>366<\/th>\n      <th>367<\/th>\n      <th>368<\/th>\n      <th>369<\/th>\n      <th>370<\/th>\n      <th>371<\/th>\n      <th>372<\/th>\n      <th>373<\/th>\n      <th>374<\/th>\n      <th>375<\/th>\n      <th>376<\/th>\n      <th>377<\/th>\n      <th>378<\/th>\n      <th>379<\/th>\n      <th>380<\/th>\n      <th>381<\/th>\n      <th>382<\/th>\n      <th>383<\/th>\n      <th>384<\/th>\n      <th>385<\/th>\n      <th>386<\/th>\n      <th>387<\/th>\n      <th>388<\/th>\n      <th>389<\/th>\n      <th>390<\/th>\n      <th>391<\/th>\n      <th>392<\/th>\n      <th>393<\/th>\n      <th>394<\/th>\n      <th>395<\/th>\n      <th>396<\/th>\n      <th>397<\/th>\n      <th>398<\/th>\n      <th>399<\/th>\n      <th>400<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

#### List with final recoded genotypes (complete cases only)


```r
gt_nomiss_list = purrr::map(gt_list, function(BIN_LENGTH_DF){
  BIN_LENGTH_DF %>%
    dplyr::filter(complete.cases(.))
})
```


## Simulate phenotype

### Extract samples of genotypes


```r
# Set directory
out_dir_test = file.path(lts_dir, "association_testing/20211027_test")
```

These must be written to a file because `PhenotypeSimulator` reads delimited genotypes from files.


```r
# Get random 10 loci
set.seed(10)
n_loci = 10
# NOTE: PhenotypeSimulator::readStandardGenotypes states that the genotype file must
# delim: a [delimter]-delimited file of [(NrSNPs+1) x (NrSamples+1)] genotypes with the snpIDs in the first column and the sampleIDs in the first row and genotypes encoded as numbers between 0 and 2 representing the (posterior) mean genotype, or dosage of the minor allele.

sample_gts = purrr::map(seq_along(gt_nomiss_list), function(COUNTER){
  
  # Pull out random SNPs
  snp_sample = gt_nomiss_list[[COUNTER]] %>% 
    dplyr::slice_sample(n = n_loci) %>% 
    dplyr::arrange(CHROM, BIN_START) %>% 
    # create locus column
    dplyr::mutate(LOCUS = paste(CHROM, BIN_START, sep = ":")) %>% 
    # remove superfluous columns
    dplyr::select(-c(CHROM, BIN, BIN_START, BIN_END)) %>% 
    # recode back to 0,1,2
    dplyr::mutate(dplyr::across(-LOCUS,
                                ~dplyr::recode(.x,
                                               `-1` = 0,
                                               `0` = 1,
                                               `1` = 2))) %>% 
    # reorder columns
    dplyr::select(LOCUS, everything())
  
})
names(sample_gts) = names(gt_nomiss_list)

purrr::map(seq_along(sample_gts), function(COUNTER){
  # save to file
  bin_length = names(sample_gts)[COUNTER]
  out_file = file.path(out_dir_test, "target_loci", paste(bin_length, ".csv", sep = ""))
  readr::write_csv(sample_gts[[COUNTER]], out_file)
})

saveRDS(sample_gts, file.path(out_dir_test, "target_loci/sample_gts.rds"))
```




### Simulate phenotype


```r
sim_path = file.path(out_dir_test, "simulated_phenotypes/20211103_sim_phenos.rds")
```


```r
set.seed(5671)
# N samples
N = n_samples
# N phenotypes
P = 1
# Porportion of total genetic variance
genVar = 0.5
# Proportion of genetic variance of genetic variant effects
h2s = 1
# Proportion of total noise variance
noiseVar = 0.5
# Proportion of noise variance of observational noise effects
phi = 1

sim_phenos = purrr::map(seq_along(sample_gts), function(COUNTER){
  # get sample file path
  bin_length = names(sample_gts)[COUNTER]
  gt_sample_file = file.path(out_dir_test, "target_loci", paste(bin_length, ".csv", sep = ""))
    
  sim_pheno = PhenotypeSimulator::runSimulation(N = N, P = P, 
                                                genVar = genVar, h2s = h2s, 
                                                noiseVar = noiseVar, phi = phi,
                                                cNrSNP = n_loci,
                                                genotypefile = gt_sample_file,
                                                format = "delim",
                                                genoDelimiter = ",")
  
  return(sim_pheno)
})
names(sim_phenos) = names(sample_gts)

saveRDS(sim_phenos, sim_path)
```



## Reformat genotypes for GridLMM

### Complete cases


```r
final_nomiss = purrr::map(seq_along(gt_nomiss_list), function(COUNTER){
  out = list()
  
  # Genotypes
  out[["genotypes"]] = gt_nomiss_list[[COUNTER]] %>% 
    dplyr::select(-c(CHROM, BIN, BIN_START, BIN_END)) %>% 
    # convert to matrix
    as.matrix(.) %>% 
    # transpose to put samples as rows
    t(.) %>% 
    # convert to data frame
    as.data.frame(.)
  
  # Positions
  out[["positions"]] = gt_nomiss_list[[COUNTER]] %>% 
    dplyr::select(CHROM, BIN_START, BIN_END)
  
  # Phenotypes
  out[["phenotypes"]] = data.frame(SAMPLE = rownames(sim_phenos[[COUNTER]]$phenoComponentsFinal$Y),
                                   Y = sim_phenos[[COUNTER]]$phenoComponentsFinal$Y %>%
                                     as.numeric()
                                   )
    
  return(out)
})
names(final_nomiss) = names(gt_nomiss_list)
```

### With NAs


```r
final_wimiss = purrr::map(seq_along(gt_list), function(COUNTER){
  out = list()
  
  # Genotypes
  out[["genotypes"]] = gt_list[[COUNTER]] %>% 
    dplyr::select(-c(CHROM, BIN, BIN_START, BIN_END)) %>% 
    # convert to matrix
    as.matrix(.) %>% 
    # transpose to put samples as rows
    t(.) %>% 
    # convert to data frame
    as.data.frame(.)
  
  # Positions
  out[["positions"]] = gt_list[[COUNTER]] %>% 
    dplyr::select(CHROM, BIN_START, BIN_END)
  
  # Phenotypes
  out[["phenotypes"]] = data.frame(SAMPLE = rownames(sim_phenos[[COUNTER]]$phenoComponentsFinal$Y),
                                   Y = sim_phenos[[COUNTER]]$phenoComponentsFinal$Y %>%
                                     as.numeric()
                                   )
    
  return(out)
})
names(final_wimiss) = names(gt_list)
```

## Test GridLMM

### No missing

#### Run GWAS


```r
test_out_nomiss = file.path(out_dir_test, "gwas_results", "gwas_results_nomiss.rds")
```


```r
gwas_tests_nomiss = purrr::map(final_nomiss, function(BIN_LENGTH){
  run_gwas(d = BIN_LENGTH$genotypes,
           m = BIN_LENGTH$positions,
           p = BIN_LENGTH$phenotypes)
})

saveRDS(gwas_tests_nomiss, test_out_nomiss)
```



#### Plot


```r
# Create custom Manhattan plot

purrr::map(seq_along(gwas_tests_nomiss), function(COUNTER){
  # Get bin length
  BIN_LENGTH = names(gwas_tests_nomiss)[COUNTER] %>% 
    as.numeric(.)
  
  # Clean data frame
  test_results = gwas_tests_nomiss[[COUNTER]]$results %>% 
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
    dplyr::mutate(TARGET = dplyr::if_else(LOCUS %in% sample_gts[[COUNTER]]$LOCUS,
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
    ylab("-log10(p-value)") + 
    geom_hline(yintercept = significance_line, colour = "#1effbc", linetype = "dashed") +
    geom_hline(yintercept = suggestive_line, colour = "#5c95ff", linetype = "dashed")
  
  ggplotly(p1, tooltip = c("CHROM", "BIN_START", "BIN_END"))

})
#> [[1]]
#> 
#> [[2]]
#> 
#> [[3]]
#> 
#> [[4]]
```

### Include missing genotypes

#### Run GWAS


```r
test_out_wimiss = file.path(out_dir_test, "gwas_results", "gwas_results_wimiss.rds")
```


```r
gwas_tests_wimiss = purrr::map(final_wimiss, function(BIN_LENGTH){
  run_gwas(d = BIN_LENGTH$genotypes,
           m = BIN_LENGTH$positions,
           p = BIN_LENGTH$phenotypes)
})

saveRDS(gwas_tests_wimiss, test_out_wimiss)
```



#### Plot


```r
# Create custom Manhattan plot
gwas_pal = c("#2B2D42", "#F7B267", "#F25C54")
names(gwas_pal) = c("target", "even chr", "odd chr")
significance_line = 3.6
suggestive_line = 2.9

purrr::map(seq_along(gwas_tests_wimiss), function(COUNTER){
  # Get bin length
  BIN_LENGTH = names(gwas_tests_wimiss)[COUNTER] %>% 
    as.numeric(.)
  
  # Clean data frame
  test_results = gwas_tests_wimiss[[COUNTER]]$results %>% 
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
    dplyr::mutate(TARGET = dplyr::if_else(LOCUS %in% sample_gts[[COUNTER]]$LOCUS,
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
    ylab("-log10(p-value)") + 
    geom_hline(yintercept = significance_line, colour = "#1effbc", linetype = "dashed") +
    geom_hline(yintercept = suggestive_line, colour = "#5c95ff", linetype = "dashed")
  
  ggplotly(p1, tooltip = c("CHROM", "BIN_START", "BIN_END"))

})
#> [[1]]
#> 
#> [[2]]
#> 
#> [[3]]
#> 
#> [[4]]
```

## Phenotype data

First 400 phenotype data: <https://github.com/brettellebi/somites/tree/master/data/20210917_First400_F2_DF.xlsx>.

>I also attach here the DataFrame of the phenotyping of the first 400 F2s (First400_F2_DF.xlsx)...

>Of interest for us for the association testing are the:
 intercept period -> intercept in the table
 mean period -> mean in the table
 
Snakemake rules for running the GWAS over these phenotypes: <https://github.com/brettellebi/somites/blob/master/workflow/rules/07_assocation_testing.smk>
 
### Results

#### Read in files


```r
target_phenos = c("mean", "intercept")
names(target_phenos) = target_phenos
results_dir = "/nfs/research/birney/users/ian/somites/association_testing/20211104_true/results"

gwas_true = purrr::map(target_phenos, function(PHENO){
  purrr::map(bin_lengths, function(BIN_LENGTH){
    readRDS(file.path(results_dir, PHENO, paste(BIN_LENGTH, ".rds", sep = "")))
  })
})
```

#### Manhattans


```r
# Process data
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

# Plot
purrr::map(seq_along(plot_dat), function(COUNTER_1){
  purrr::map(seq_along(plot_dat[[COUNTER_1]]), function(COUNTER_2){
    # Get bin length
    BIN_LENGTH = names(plot_dat[[COUNTER_1]])[COUNTER_2] %>% 
      as.numeric()
    # Plot
    plot_int_man(plot_dat[[COUNTER_1]][[COUNTER_2]],
                 phenotype = names(plot_dat)[COUNTER_1],
                 bin_length = names(plot_dat[[COUNTER_1]])[COUNTER_2],
                 gwas_pal = gwas_pal[2:3],
                 med_chr_lens = med_chr_lens)
    
  })
})
#> [[1]]
#> [[1]][[1]]
#> 
#> [[1]][[2]]
#> 
#> [[1]][[3]]
#> 
#> [[1]][[4]]
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 
#> [[2]][[2]]
#> 
#> [[2]][[3]]
#> 
#> [[2]][[4]]
```


