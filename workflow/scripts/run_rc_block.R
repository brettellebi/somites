#!/usr/bin/env Rscript

library(ViteRbi)
source(snakemake@input[["source_code"]])

dp_files = snakemake@input[["dp_files"]]

pin = "all"
bin1 =collect_chunked_data(dp_files, 1, length(dp_files))
colnames(bin1)[1] = "sample"
write.table(bin1, file=snakemake@output[[1]], sep="\t", row.names=F, quote=F)

