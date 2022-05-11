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
#GENOS = "/hps/nobackup/birney/users/ian/somites/hmm_out/F2/hdrr/hmmlearn_true/None/5000/0.8.csv"
#PHENOS = "/hps/nobackup/birney/users/ian/somites/phens/inv_norm/intercept.phen"
#COV = 0.8
#MAX_READS = "None"

## True
GENOS = snakemake@input[["genos"]]
PHENOS = snakemake@input[["phenos"]]
COV = snakemake@params[["cov"]]
MAX_READS = snakemake@params[["max_reads"]]
BIN_LENGTH = as.numeric(snakemake@params[["bin_length"]])
PROP_SITES_TOT_TOP = snakemake@output[["prop_sites_total_top"]]
PROP_SITES_TOT_BOTTOM = snakemake@output[["prop_sites_total_bottom"]]
KARYOPLOT_NOMISS_TOP = snakemake@output[["karyoplot_no_missing_top"]]
KARYOPLOT_NOMISS_BOTTOM = snakemake@output[["karyoplot_no_missing_bottom"]]
KARYOPLOT_WIMISS_TOP = snakemake@output[["karyoplot_with_missing_top"]]
KARYOPLOT_WIMISS_BOTTOM = snakemake@output[["karyoplot_with_missing_bottom"]]


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

# Genos

df = readr::read_csv(GENOS,
                     col_types = "iiiiiidi") %>% 
  # add key variables
  dplyr::mutate(BIN_START = (BIN * BIN_LENGTH) + 1,
                BIN_END = ((BIN + 1) * BIN_LENGTH)) 

# Phenos

phenos = readr::read_tsv(PHENOS,
                         col_names = c("FID", "IID", "intercept")) %>% 
  dplyr::arrange(dplyr::desc(intercept)) %>% 
  tidyr::drop_na()

phen_top = phenos %>% 
  dplyr::slice_head(n = 50) %>% 
  dplyr::pull(IID)

phen_bottom = phenos %>% 
  dplyr::slice_tail(n = 50) %>% 
  dplyr::pull(IID)

# Filter df

df_top = df %>% 
  dplyr::filter(SAMPLE %in% phen_top) %>% 
  dplyr::arrange(match(SAMPLE, phen_top)) %>% 
  # make factor to order
  dplyr::mutate(SAMPLE = factor(SAMPLE, levels = phen_top))

df_bottom = df %>% 
  dplyr::filter(SAMPLE %in% phen_bottom) %>% 
  dplyr::arrange(match(SAMPLE, phen_bottom)) %>% 
  # make factor to order
  dplyr::mutate(SAMPLE = factor(SAMPLE, levels = phen_bottom))

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
# Total sites covered by each state top
#######################

# All sites

prop_sites_tot = df_top %>% 
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

ggsave(PROP_SITES_TOT_TOP,
       prop_sites_tot,
       device = "png",
       width = 10,
       height = 5.8,
       units = "in",
       dpi = 400)

#######################
# Total sites covered by each state bottom
#######################

# All sites

prop_sites_tot = df_bottom %>% 
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

ggsave(PROP_SITES_TOT_BOTTOM,
       prop_sites_tot,
       device = "png",
       width = 10,
       height = 5.8,
       units = "in",
       dpi = 400)

#######################
# Karyoplot no missing top
#######################

# Create custom genome 

med_genome = regioneR::toGRanges(med_chr_lens)

# Convert data to list of block boundaries for each LANE

block_bounds_list = df_top %>% 
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

png(file=KARYOPLOT_NOMISS_TOP,
    width=7800,
    height=3000,
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Top 50 samples for split inverse-normalised period intercept",
                                  "\nEmission (co)variances: ", 
                                  COV,
                                  "\nMax reads per bin: ",
                                  MAX_READS,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=0.5)

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

#######################
# Karyoplot no missing bottom
#######################

# Create custom genome 

med_genome = regioneR::toGRanges(med_chr_lens)

# Convert data to list of block boundaries for each LANE

block_bounds_list = df_bottom %>% 
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

png(file=KARYOPLOT_NOMISS_BOTTOM,
    width=7800,
    height=3000,
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Bottom 50 samples for split inverse-normalised period intercept",
                                  "\nEmission (co)variances: ", 
                                  COV,
                                  "\nMax reads per bin: ",
                                  MAX_READS,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=0.5)

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


#######################
# Karyoplot WITH missing top
#######################

## Process
block_bounds_list = df_top %>% 
  # loop over LANE
  split(., f = .$SAMPLE) %>% 
  purrr::map(., function(SAMPLE){
    
    STRAIN = unique(SAMPLE$SAMPLE)
    # Create list of possible bins
    poss_bins = purrr::map(med_chr_lens$chr, function(CHROM){
      # Get chr end
      CHR_END = med_chr_lens %>% 
        dplyr::filter(chr == CHROM) %>% 
        dplyr::pull(end) %>% 
        as.numeric()
      # Get bin starts
      out = tibble::tibble(CHROM = as.numeric(CHROM),
                           BIN_START = seq(from = 1, to = CHR_END, by = BIN_LENGTH),
                           BIN_END = BIN_START + BIN_LENGTH - 1
      )
      # Adjust final bin end 
      out[nrow(out), "BIN_END"] = CHR_END
      
      return(out)
    }) %>% 
      dplyr::bind_rows()
    
    
    # Bind DF
    new_df = dplyr::left_join(poss_bins,
                              SAMPLE %>% 
                                dplyr::select(CHROM, BIN_START, BIN_END, STATE),
                              by = c("CHROM", "BIN_START", "BIN_END")) %>% 
      # replace NAs with `UNCLASSIFIED`
      dplyr::mutate(STATE = as.character(STATE),
                    STATE = STATE %>% 
                      tidyr::replace_na("UNCLASSIFIED"),
                    # add STRAIN
                    SAMPLE = STRAIN) %>% 
      # add COLOUR
      dplyr::mutate(COLOUR = dplyr::recode(STATE,
                                           !!!pal_hom_het_2))
    
    
  })

# Plot karyoplot

png(file=KARYOPLOT_WIMISS_TOP,
    width=7800,
    height=3000,
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Top 50 samples for split inverse-normalised period intercept",
                                  "Emission (co)variances: ", 
                                  COV,
                                  "\nMax reads per bin: ",
                                  MAX_READS,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=0.5)


# Add rectangles in loop
counter = 0
purrr::map(block_bounds_list, function(SAMPLE){
  # Add to counter
  counter <<- counter + 1
  # Add rectangles
  karyoploteR::kpRect(kp,
                      chr = SAMPLE$CHROM,
                      x0 = SAMPLE$BIN_START,
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

#######################
# Karyoplot WITH missing bottom
#######################

## Process
block_bounds_list = df_bottom %>% 
  # loop over LANE
  split(., f = .$SAMPLE) %>% 
  purrr::map(., function(SAMPLE){
    
    STRAIN = unique(SAMPLE$SAMPLE)
    # Create list of possible bins
    poss_bins = purrr::map(med_chr_lens$chr, function(CHROM){
      # Get chr end
      CHR_END = med_chr_lens %>% 
        dplyr::filter(chr == CHROM) %>% 
        dplyr::pull(end) %>% 
        as.numeric()
      # Get bin starts
      out = tibble::tibble(CHROM = as.numeric(CHROM),
                           BIN_START = seq(from = 1, to = CHR_END, by = BIN_LENGTH),
                           BIN_END = BIN_START + BIN_LENGTH - 1
      )
      # Adjust final bin end 
      out[nrow(out), "BIN_END"] = CHR_END
      
      return(out)
    }) %>% 
      dplyr::bind_rows()
    
    
    # Bind DF
    new_df = dplyr::left_join(poss_bins,
                              SAMPLE %>% 
                                dplyr::select(CHROM, BIN_START, BIN_END, STATE),
                              by = c("CHROM", "BIN_START", "BIN_END")) %>% 
      # replace NAs with `UNCLASSIFIED`
      dplyr::mutate(STATE = as.character(STATE),
                    STATE = STATE %>% 
                      tidyr::replace_na("UNCLASSIFIED"),
                    # add STRAIN
                    SAMPLE = STRAIN) %>% 
      # add COLOUR
      dplyr::mutate(COLOUR = dplyr::recode(STATE,
                                           !!!pal_hom_het_2))
    
    
  })

# Plot karyoplot

png(file=KARYOPLOT_WIMISS_BOTTOM,
    width=7800,
    height=3000,
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Bottom 50 samples for split inverse-normalised period intercept",
                                  "Emission (co)variances: ", 
                                  COV,
                                  "\nMax reads per bin: ",
                                  MAX_READS,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=0.5)


# Add rectangles in loop
counter = 0
purrr::map(block_bounds_list, function(SAMPLE){
  # Add to counter
  counter <<- counter + 1
  # Add rectangles
  karyoploteR::kpRect(kp,
                      chr = SAMPLE$CHROM,
                      x0 = SAMPLE$BIN_START,
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
