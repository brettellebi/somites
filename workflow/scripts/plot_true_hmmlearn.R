# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(karyoploteR)
library(ggbeeswarm)

# Get variables

## Debug
#BIN_LENGTH = as.numeric("5000")
#IN_FILE = "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn_true/None/5000/C.csv"
#COV = "C"
#MAX_READS = "None"

## True
IN_FILE = snakemake@input[[1]]
COV = snakemake@params[["cov"]]
MAX_READS = snakemake@params[["max_reads"]]
BIN_LENGTH = as.numeric(snakemake@params[["bin_length"]])
SCATTER = snakemake@output[["scatter"]]
BASE_COV_TOT = snakemake@output[["base_cov_total"]]
PROP_SITES_TOT = snakemake@output[["prop_sites_total"]]
KARYOPLOT_NOMISS = snakemake@output[["karyoplot_no_missing"]]
KARYOPLOT_WIMISS = snakemake@output[["karyoplot_with_missing"]]


######################
# Number of states
######################

N_STATES = 3

######################
# Palette
######################

pal_hom_het_2 = c("#43AA8B", "#000022", "#DE3C4B", "#FBF5F3")
names(pal_hom_het_2) = c(0:2, "UNCLASSIFIED")

pal_hom_het_2_lines = c(karyoploteR::darker(pal_hom_het_2[1], 100),
                        karyoploteR::lighter(pal_hom_het_2[2], 100),
                        karyoploteR::darker(pal_hom_het_2[3], 100),
                        karyoploteR::darker(pal_hom_het_2[4]))
names(pal_hom_het_2_lines) = c(0:2, "UNCLASSIFIED")


# Recode vector

recode_vec = c(`0` = "Homozygous Cab",
               `1` = "Heterozygous",
               `2` = "Homozygous Kaga")


# Set states to loop over

states = 0:(N_STATES - 1)

######################
# Read in data
######################

df = readr::read_csv(IN_FILE,
                     col_types = "iiiiiidi") %>% 
  # add key variables
  dplyr::mutate(BIN_START = (BIN * BIN_LENGTH) + 1,
                BIN_END = ((BIN + 1) * BIN_LENGTH)) 

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

#######################
# Total sites covered by each state
#######################

# All sites

prop_sites_tot = df %>% 
    # get counts of sites per LANE and state
    dplyr::group_by(SAMPLE, STATE) %>% 
    dplyr::count() %>% 
    # spread to one row per LANE
    tidyr::pivot_wider(id_cols = SAMPLE, names_from = STATE, values_from = n) %>% 
    # calculate frequencies of states per LANE
    dplyr::mutate(TOTAL = sum(`0`, `1`, `2`),
                  FREQ_0 = `0` / TOTAL,
                  FREQ_1 = `1` / TOTAL,
                  FREQ_2 = `2` / TOTAL) %>% 
    # gather
    tidyr::pivot_longer(cols = starts_with("FREQ_"),
                        names_to = "STATE",
                        names_prefix = "FREQ_",
                        values_to = "FREQ") %>% 
    # order STATE and recode with meaning
    dplyr::mutate(STATE = factor(STATE, levels = states),
                  STATE_RECODE = dplyr::recode(STATE, !!!recode_vec)) %>% 
    # plot
    ggplot(aes(STATE_RECODE, FREQ, colour = STATE, fill = STATE)) +
    geom_violin() +
    geom_boxplot(width = .5) +
    ggbeeswarm::geom_quasirandom(color="#7D8491", size=0.4, alpha=0.9) +
    theme_bw() +
    scale_colour_manual(values = pal_hom_het_2_lines) +
    scale_fill_manual(values = pal_hom_het_2) +
    guides(colour = "none", fill = "none") +
    xlab("Genotype") +
    ylab("Frequency") +
    ggtitle(paste("Emission (co)variances: ", 
                  COV,
                  "\nMax reads per bin: ",
                  MAX_READS,
                  "\nBin length: ",
                  BIN_LENGTH,
                  sep = ""))

ggsave(PROP_SITES_TOT,
       prop_sites_tot,
       device = "png",
       width = 10,
       height = 5.8,
       units = "in",
       dpi = 400)

#######################
# Karyoplots
#######################

# Create custom genome 

med_genome = regioneR::toGRanges(med_chr_lens)

# Convert data to list of block boundaries for each LANE

block_bounds_list = df %>% 
  # loop over LANE
  split(., f = .$SAMPLE) %>% 
  purrr::map(., function(SAMPLE){
    # loop over CHR
    SAMPLE %>% 
      split(., f = .$CHROM) %>% 
      purrr::map(., function(CHROM){
        # Get lengths of each contiguous state
        cont_len = rle(CHROM$STATE)
        
        # Get cumulative sum of those lengths
        cum_blocks = cumsum(cont_len$lengths)
        
        # Get rows that correspond to block changes
        block_bounds = CHROM[cum_blocks, ] %>% 
          # Add end of previous block
          dplyr::mutate(END_PREV = dplyr::lag(BIN_END)) %>% 
          # Replace the NA in the first row with `1`
          dplyr::mutate(END_PREV = tidyr::replace_na(END_PREV, 1)) %>% 
          # Add colour
          dplyr::mutate(COLOUR = dplyr::recode(STATE,
                                               !!!pal_hom_het_2[-which(names(pal_hom_het_2) == "UNCLASSIFIED")])) 
        
      }) %>% 
      dplyr::bind_rows()
    
  })

# Extract y cutoff points for each lane

lane_cutoffs = cut(0:1, breaks = length(block_bounds_list), dig.lab = 7) %>% 
  levels(.) %>% 
  data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", .) ),
             upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", .) )) %>% 
  dplyr::arrange(dplyr::desc(lower))

# Plot karyoplot WITH NO missing blocks

png(file=KARYOPLOT_NOMISS,
    width=7800,
    height=23400,
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Emission (co)variances: ", 
                                  COV,
                                  "\nMax reads per bin: ",
                                  MAX_READS,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=4)

# Add rectangles in loop
counter = 0
purrr::map(block_bounds_list, function(SAMPLE){
  # Add to counter_B
  counter <<- counter + 1
  # Add rectangles
  karyoploteR::kpRect(kp,
                      chr = SAMPLE$CHROM,
                      x0 = SAMPLE$END_PREV,
                      x1 = SAMPLE$BIN_END,
                      y0 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(lower),
                      y1 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(upper),
                      col = SAMPLE$COLOUR,
                      border = NA)
  # Add axis label
  karyoploteR::kpAddLabels(kp, labels = unique(SAMPLE$SAMPLE),
                           r0 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(lower),
                           r1 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(upper),
                           cex = 0.5)
})


dev.off()  
