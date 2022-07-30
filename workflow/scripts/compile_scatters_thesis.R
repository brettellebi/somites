# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Set variables

## Debug
IN = list(
  "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn/None/5000/A.csv",
  "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn/None/5000/B.csv",
  "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn/None/5000/D.csv",
  "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn/None/5000/C.csv",
  "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn/None/5000/F.csv",
  "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn/None/5000/G.csv"
)
OUT_PNG = here::here("book/plots/thesis/scatter_collage.png")
OUT_PDF = here::here("book/plots/thesis/scatter_collage.pdf")

## True

IN = snakemake@input
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]

###################
# Extra variables
###################

# Plot only chr 18 as an example
TARGET_CHROM = 18

BIN_LENGTH = 5000

states = 0:2

pal_hom_het_2 = c("#43AA8B", "#000022", "#DE3C4B")
names(pal_hom_het_2) = states

# Set recode vector

recode_vec = c(`0` = "Homozygous Cab",
               `1` = "Heterozygous",
               `2` = "Homozygous Kaga")

######################
# Function for titles
######################

get_titles = function(MOD){
  if (MOD == "A"){
    mod = "Standard HMM"
  } else if (MOD == "B"){
    mod = "Fixed transition probabilities"
  } else if (MOD == "C"){
    mod = "Fixed transition probabilities\nEmission variances of 0.33"
  } else if (MOD == "D"){
    mod = "Fixed transition probabilities\nEmission variances of 0.01"
  } else if (MOD == "E"){
    mod = "Error states (0.15 prob of entering error state)\nVariances: 0.2 for true state, 1 for error state"
  } else if (MOD == "F"){
    mod = "Fixed transition probabilities\nEmission variances of 0.8"
  } else if (MOD == "G"){
    mod = "Fixed transition probabilities\nEmission variances of 1"
  }
  return(mod)
}


###################
# Read in data
###################

dat_list = purrr::map(IN, function(DF){
  DF %>% 
    readr::read_csv(., col_types = "iiiiiidi") %>% 
    # add key variables
    dplyr::mutate(BIN_START = (BIN * BIN_LENGTH) + 1,
                  BIN_END = ((BIN + 1) * BIN_LENGTH)) 
})

# Set names
names(dat_list) = unlist(IN) %>% 
  basename() %>% 
  stringr::str_remove(".csv")

# Read in total medaka genome count

## Get chromosome lengths
med_chr_lens = read.table("data/Oryzias_latipes.ASM223467v1.dna.toplevel.fa_chr_counts.txt",
                          col.names = c("chr", "end"))
## Add start
med_chr_lens$start = 1
## Reorder
med_chr_lens = med_chr_lens %>% 
  dplyr::select(chr, start, end) %>% 
  # remove MT
  dplyr::filter(chr != "MT")
## Total HdrR sequence length
total_hdrr_bases = sum(med_chr_lens$end)


###################
# Plot
###################

# Make list of plots

counter = 0
fig_list = purrr::map(dat_list, function(DF){
  counter <<- counter + 1
  # Get title
  mod = get_titles(names(dat_list)[counter])
  # Plot
  scatter = DF %>%
    # get Mb
    dplyr::mutate(BIN_MB = BIN_START/1e6) %>% 
    dplyr::filter(CHROM == TARGET_CHROM) %>% 
    dplyr::mutate(STATE = factor(STATE, levels = states)) %>% 
    ggplot() +
    geom_point(aes(BIN_MB, PROP_KAGA, colour = STATE),
               size = 0.5) +
    scale_colour_manual(values = pal_hom_het_2) +
    facet_grid(rows = vars(SAMPLE)) +
    cowplot::theme_cowplot(rel_large = 13/14,
                           rel_small = 11/14) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    ggtitle(mod) +
    xlab("bin start (Mb)") + 
    ylab("frequency of SNPs within bin\nthat align with the Kaga allele")
})

# Compile

legend = cowplot::get_legend(fig_list[[1]] + 
                               guides(color = guide_legend(title = "HMM state")))

out = cowplot::plot_grid(fig_list[[1]] +
                            guides(colour = "none") +
                            theme(legend.position='none',
                                  axis.title.y = element_blank()),
                          fig_list[[2]] +
                           guides(colour = "none") +
                           theme(legend.position='none') +
                            theme(axis.title.y = element_blank()),
                          fig_list[[3]]  +
                           guides(colour = "none") +
                           theme(legend.position='none'),
                          fig_list[[4]]  +
                           guides(colour = "none") +
                           theme(legend.position='none') +
                            theme(axis.title.y = element_blank()),
                          fig_list[[5]]  +
                           guides(colour = "none") +
                           theme(legend.position='none',
                                 axis.title.y = element_blank()),
                          fig_list[[6]]  +
                           guides(colour = "none") +
                           theme(legend.position='none') +
                            theme(axis.title.y = element_blank()),
                          align = "hv", ncol = 2, axis = "tblr",
                         labels = c("A", "B", "C", "D", "E", "F"))

# Add legend
out = cowplot::plot_grid(out,
                         cowplot::plot_grid(NULL, legend, NULL, ncol = 1),
                         rel_widths=c(1, 0.2))
    

ggsave(OUT_PNG,
       out,
       device = "png",
       width = 12,
       height = 15,
       units = "in",
       dpi = 400)

ggsave(OUT_PDF,
       out,
       device = "pdf",
       width = 12,
       height = 15,
       units = "in",
       dpi = 400)
