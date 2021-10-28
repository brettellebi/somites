# Association testing with `GridLMM`

Example files for GridLMM:
`/nfs/research/birney/users/ian/somites/gridLMM_example`, copied from:
`/nfs/research/birney/users/tomas/scratch/gridlmm_inputs`

Should be fairly self explanatory, but:

* `all_genotypes.txt`: coded -1,0,1 (AA, AB, BB), rows are samples, columns are "SNPs" (recombination blocks).

* `all_positions.txt`: genome coordinates (chr, start, end) **NB**. In the example I simply put an index which related to a different file (i.e. the merged recombination blocks) but this could (should) have had the real positions.

* `phenotypes.txt`: this has sample_id (or file name) and phenotype.

**NB**. Super important that the genotype file and phenotype file are ordered in the same order, i.e. row1 in genotype file in the same sample as row one is phenotype file and so on.

Btw: to successfully install GridLMM i needed to not install the vignettes, I used:
```
devtools::install_github('deruncie/GridLMM', build_opts = c("--no-resave-data", "--no-manual"),force = TRUE,build_vignettes = FALSE)
```