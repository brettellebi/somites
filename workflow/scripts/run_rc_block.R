

library(ViteRbi)
source("rc_block_hmm.R")

dp_file_path = "../dpAB_all/"
dp_files = paste(dp_file_path, dir(dp_file_path), sep="")
dp_files = as.list(dp_files)

pin = "all"
bin1 =collect_chunked_data(dp_files, 1, length(dp_files))
colnames(bin1)[1] = "sample"
write.table(bin1, file=paste(pin, "_hmm_output.txt", sep=""), sep="\t", row.names=F, quote=F)

