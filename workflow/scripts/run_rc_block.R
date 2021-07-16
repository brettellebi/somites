#!/usr/bin/env Rscript

library(ViteRbi)
source(snakemake@input[["source_code"]])

dp_files = snakemake@input[["dp_files"]]

# Set bin length
bin_len = as.integer(snakemake@params[["bin_length"]])

pin = "all"
bin1 =collect_chunked_data(dp_files, 1, length(dp_files), bin_len = bin_len)
colnames(bin1)[1] = "sample"
write.table(bin1, file=snakemake@output[[1]], sep="\t", row.names=F, quote=F)

