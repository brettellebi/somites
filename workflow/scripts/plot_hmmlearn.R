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
BIN_LENGTH = as.numeric("5000")
IN_FILE = "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn/None/5000/E.csv"
MOD = "E"
MAX_READS = "None"

## True
IN_FILE = snakemake@input[[1]]
MOD = snakemake@params[["mod"]]
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

N_STATES = c("A" = 3,
             "B" = 3,
             "C" = 3,
             "D" = 3,
             "E" = 5,
             "F" = 3,
             "G" = 3)

N_STATES = N_STATES[MOD]

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

# Extend palette for "error" states

if (N_STATES == 5){
  pal_hom_het_2 = c("#43AA8B", karyoploteR::lighter("#43AA8B"),
                    "#000022",
                    karyoploteR::lighter("#DE3C4B"), "#DE3C4B", 
                    "#FBF5F3")
  pal_hom_het_2_lines = c(karyoploteR::darker(pal_hom_het_2[1], 100),
                          karyoploteR::darker(pal_hom_het_2[2], 100),
                          karyoploteR::lighter(pal_hom_het_2[3], 100),
                          karyoploteR::darker(pal_hom_het_2[4], 100),
                          karyoploteR::darker(pal_hom_het_2[5], 100))
  names(pal_hom_het_2) = c(0:4, "UNCLASSIFIED")
}

# Recode vector

if (N_STATES == 5){
  recode_vec = c(`0` = "Homozygous Cab",
                 `1` = "Cab error",
                 `2` = "Heterozygous",
                 `3` = "Kaga error",
                 `4` = "Homozygous Kaga")
} else {
  recode_vec = c(`0` = "Homozygous Cab",
                 `1` = "Heterozygous",
                 `2` = "Homozygous Kaga")
}


######################
# Elaborate on mods
######################

if (MOD == "A"){
  mod = "Standard HMM"
} else if (MOD == "B"){
  mod = "Fixed transition probabilities"
} else if (MOD == "C"){
  mod = "Fixed transition probabilities and wide, even covariances (0.33)"
} else if (MOD == "D"){
  mod = "Fixed transition probabilities and narrow, even covariances (0.01)"
} else if (MOD == "E"){
  mod = "Error states (0.15 prob of entering error state)\nVariances: 0.2 for true state, 1 for error state"
} else if (MOD == "F"){
  mod = "Fixed transition probabilities and wide, even covariances (0.8)"
} else if (MOD == "G"){
  mod = "Fixed transition probabilities and wide, even covariances (1)"
}

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


######################
# Scatter plot
######################

TARGET_CHROM = 18
scatter = df %>%
  dplyr::filter(CHROM == TARGET_CHROM) %>% 
  dplyr::mutate(STATE = factor(STATE, levels = states)) %>% 
  ggplot() +
  geom_point(aes(BIN_START, PROP_KAGA, colour = STATE),
             size = 0.5) +
  scale_colour_manual(values = pal_hom_het_2) +
  facet_grid(rows = vars(SAMPLE)) +
  theme_bw() +
  ggtitle(paste("Modification: ",
                mod,
                "\nChromosome: ",
                TARGET_CHROM))

ggsave(SCATTER,
       scatter,
       device = "png",
       width = 10,
       height = 5.8,
       units = "in",
       dpi = 400)

########################
## Total number of bases covered by each state
########################
#

#
#names(states) = states
#
## All sites
#
### Process
#base_cov_df = df %>% 
#  split(., f = .$SAMPLE) %>% 
#  purrr::map(., function(SAMPLE){
#    # convert to ranges object
#    lane_ranges = GenomicRanges::makeGRangesFromDataFrame(SAMPLE,
#                                                          keep.extra.columns = T,
#                                                          ignore.strand = T,
#                                                          seqnames.field = "CHROM", 
#                                                          start.field = "BIN_START",
#                                                          end.field = "BIN_END")
#    # get total bases covered by each state
#    purrr::map_dfr(states, function(STATE){
#      lane_ranges[lane_ranges$STATE == STATE] %>% 
#        # merge contiguous ranges
#        GenomicRanges::reduce(.) %>% 
#        # get width of ranges
#        width(.) %>% 
#        # get total bases covered
#        sum(.) %>% 
#        # coerce into data frame
#        data.frame("BASES_COVERED" = .)
#    }, .id = "STATE") %>% 
#      # add FREQ column
#      dplyr::mutate(FREQ = BASES_COVERED / total_hdrr_bases) %>% 
#      # add UNCLASSIFIED row
#      tibble::add_row(STATE = "UNCLASSIFIED", 
#                      BASES_COVERED = total_hdrr_bases - sum(.$BASES_COVERED),
#                      FREQ = (total_hdrr_bases - sum(.$BASES_COVERED)) / total_hdrr_bases)
#  }
#  ) %>% 
#  dplyr::bind_rows(.id = "LANE")
#
#
### Plot
#base_cov_tot = base_cov_df %>% 
#  dplyr::mutate(STATE = factor(STATE, levels = c(states, "UNCLASSIFIED")),
#                STATE_RECODE = dplyr::recode(STATE, !!!recode_vec)) %>% 
#  # plot
#  ggplot(aes(STATE_RECODE, FREQ, colour = STATE, fill = STATE)) +
#  geom_violin() +
#  geom_boxplot(width = .3) +
#  ggbeeswarm::geom_quasirandom(color="#7D8491", size=0.4, alpha=0.9) +
#  theme_bw() +
#  scale_colour_manual(values = pal_hom_het_2_lines) +
#  scale_fill_manual(values = pal_hom_het_2) +
#  guides(colour = "none", fill = "none") +
#  xlab("Genotype") +
#  ylab("Proportion of reference bases covered") +
#  ggtitle(paste("Parameter modification: ", 
#                mod,
#                "\nMax reads per bin: ",
#                MAX_READS,
#                "\nBin length: ",
#                BIN_LENGTH,
#                sep = ""))
#
### Save
#ggsave(BASE_COV_TOT,
#       base_cov_tot,
#       device = "png",
#       width = 10,
#       height = 5.8,
#       units = "in",
#       dpi = 400)
#
#######################
# Total sites covered by each state
#######################

# All sites

if (N_STATES == 5){
  states_with_counts = df %>% 
    count(STATE) %>% 
    dplyr::pull(STATE) %>% 
    as.character()
  
  prop_sites_tot = df %>% 
    # get counts of sites per LANE and state
    dplyr::group_by(SAMPLE, STATE) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    # spread to one row per LANE
    tidyr::pivot_wider(id_cols = SAMPLE, names_from = STATE, values_from = n) %>% 
    # calculate frequencies of states per LANE
    dplyr::mutate(TOTAL = rowSums(across(all_of(states_with_counts))),
                  dplyr::across(dplyr::all_of(states_with_counts),
                                ~ .x / TOTAL,
                                .names = "FREQ_{.col}")) %>% 
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
    ggtitle(paste("Parameter modification: ", 
                  mod,
                  "\nMax reads per bin: ",
                  MAX_READS,
                  "\nBin length: ",
                  BIN_LENGTH,
                  sep = ""))
} else {
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
    ggtitle(paste("Parameter modification: ", 
                  mod,
                  "\nMax reads per bin: ",
                  MAX_READS,
                  "\nBin length: ",
                  BIN_LENGTH,
                  sep = ""))
  
}

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
    height=1500,
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Parameter modification: ", 
                                  mod,
                                  "\nMax reads per bin: ",
                                  MAX_READS,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=1/3)

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

#
## Plot karyoplot WITH missing blocks
#
### Process
#block_bounds_list = df %>% 
#  # loop over LANE
#  split(., f = .$SAMPLE) %>% 
#  purrr::map(., function(SAMPLE){
#  
#    STRAIN = unique(SAMPLE$SAMPLE)
#    # Create list of possible bins
#    poss_bins = purrr::map(med_chr_lens$chr, function(CHROM){
#      # Get chr end
#      CHR_END = med_chr_lens %>% 
#        dplyr::filter(chr == CHROM) %>% 
#        dplyr::pull(end) %>% 
#        as.numeric()
#      # Get bin starts
#      out = tibble::tibble(CHROM = as.numeric(CHROM),
#                           BIN_START = seq(from = 1, to = CHR_END, by = BIN_LENGTH),
#                           BIN_END = BIN_START + BIN_LENGTH - 1
#      )
#      # Adjust final bin end 
#      out[nrow(out), "BIN_END"] = CHR_END
#      
#      return(out)
#    }) %>% 
#      dplyr::bind_rows()
#  
#    
#    # Bind DF
#    new_df = dplyr::left_join(poss_bins,
#                              SAMPLE %>% 
#                                dplyr::select(CHROM, BIN_START, BIN_END, STATE),
#                              by = c("CHROM", "BIN_START", "BIN_END")) %>% 
#      # replace NAs with `UNCLASSIFIED`
#      dplyr::mutate(STATE = as.character(STATE),
#                    STATE = STATE %>% 
#                      tidyr::replace_na("UNCLASSIFIED"),
#                    # add STRAIN
#                    SAMPLE = STRAIN) %>% 
#      # add COLOUR
#      dplyr::mutate(COLOUR = dplyr::recode(STATE,
#                                           !!!pal_hom_het_2))
#  
#            
#  })
#
## Plot karyoplot
#
#png(file=KARYOPLOT_WIMISS,
#    width=7800,
#    height=23400,
#    units = "px",
#    res = 400)
#
## Plot ideogram
#kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)
#
## Plot title
#karyoploteR::kpAddMainTitle(kp,
#                            paste("Parameter modification: ", 
#                                  mod,
#                                  "\nMax reads per bin: ",
#                                  MAX_READS,
#                                  "\nBin length: ",
#                                  BIN_LENGTH,
#                                  sep = ""),
#                            cex=4)
#
## Add rectangles in loop
#counter = 0
#purrr::map(block_bounds_list, function(SAMPLE){
#  # Add to counter
#  counter <<- counter + 1
#  # Add rectangles
#  karyoploteR::kpRect(kp,
#                      chr = SAMPLE$CHROM,
#                      x0 = SAMPLE$BIN_START,
#                      x1 = SAMPLE$BIN_END,
#                      y0 = lane_cutoffs[counter, ] %>% 
#                        dplyr::pull(lower),
#                      y1 = lane_cutoffs[counter, ] %>% 
#                        dplyr::pull(upper),
#                      col = SAMPLE$COLOUR,
#                      border = NA)
#  # Add axis label
#  karyoploteR::kpAddLabels(kp, labels = unique(SAMPLE$SAMPLE),
#                           r0 = lane_cutoffs[counter, ] %>% 
#                             dplyr::pull(lower),
#                           r1 = lane_cutoffs[counter, ] %>% 
#                             dplyr::pull(upper),
#                           cex = 0.5)
#})
#
#
#dev.off()  
