# Send output to log

log <- file(LOG_FILE, open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

### Debug
#SITE_FILTER = "all_sites"
#BIN_LENGTH = as.numeric("20000")
#IN_FILE = file.path("/nfs/research/birney/users/ian/somites/recombination_blocks/F2",
#                    SITE_FILTER,
#                    paste(BIN_LENGTH, ".txt", sep = ""))

## True
IN_FILE = snakemake@input[[1]]
SITE_FILTER = snakemake@params[["site_filter"]]
BIN_LENGTH = as.numeric(snakemake@params[["bin_length"]])
LOG_FILE = snakemake@log[[1]]
BASE_COV_TOT = snakemake@output[["base_cov_total"]]
BASE_COV_BY_CHROM = snakemake@output[["base_cov_by_chrom"]]
PROP_SITES_TOT = snakemake@output[["prop_sites_total"]]
PROP_SITES_BY_CHROM = snakemake@output[["prop_sites_by_chrom"]]
KARYOPLOT_NOMISS = snakemake@output[["karyoplot_no_missing"]]
KARYOPLOT_WIMISS = snakemake@output[["karyoplot_with_missing"]]


#library(here)
source("book/source/03-F2_recombination.R")

# Read in data

df = readr::read_tsv(IN_FILE,
                     col_types = "ciiidii") %>% 
  # add key variables
  dplyr::mutate(LANE = basename(sample) %>% 
                  stringr::str_remove(".txt") %>% 
                  as.numeric(.),
                BIN_LENGTH = as.numeric(BIN_LENGTH),
                BIN_START = (bin - 1) * BIN_LENGTH + 1,
                BIN_END = bin * BIN_LENGTH,
                BIN_LENGTH_KB = BIN_LENGTH / 1e3,
                READS_PER_BIN = mat + pat) %>% 
  # Recode state to put in correct order]
  dplyr::mutate(state = dplyr::recode(state,
                                      `0` = 2,
                                      `1` = 1,
                                      `2` = 0))

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
# Total number of bases covered by each state
#######################

# Set states to loop over

states = 0:2
names(states) = states

# All sites

## Process
base_cov_df = df %>% 
  split(., f = .$LANE) %>% 
  purrr::map(., function(LANE){
    # convert to ranges object
    lane_ranges = GenomicRanges::makeGRangesFromDataFrame(LANE,
                                                          keep.extra.columns = T,
                                                          ignore.strand = T,
                                                          seqnames.field = "chr", 
                                                          start.field = "BIN_START",
                                                          end.field = "BIN_END")
    # get total bases covered by each state
    purrr::map_dfr(states, function(STATE){
      lane_ranges[lane_ranges$state == STATE] %>% 
        # merge contiguous ranges
        GenomicRanges::reduce(.) %>% 
        # get width of ranges
        width(.) %>% 
        # get total bases covered
        sum(.) %>% 
        # coerce into data frame
        data.frame("BASES_COVERED" = .)
    }, .id = "STATE") %>% 
      # add FREQ column
      dplyr::mutate(FREQ = BASES_COVERED / total_hdrr_bases) %>% 
      # add UNCLASSIFIED row
      tibble::add_row(STATE = "UNCLASSIFIED", 
                      BASES_COVERED = total_hdrr_bases - sum(.$BASES_COVERED),
                      FREQ = (total_hdrr_bases - sum(.$BASES_COVERED)) / total_hdrr_bases)
  }
  ) %>% 
  dplyr::bind_rows(.id = "LANE")


## Plot
base_cov_tot = base_cov_df %>% 
  dplyr::mutate(STATE = factor(STATE, levels = c(0,1,2, "UNCLASSIFIED")),
                STATE_RECODE = dplyr::recode(STATE,
                                             `0` = "Homozygous Cab",
                                             `1` = "Heterozygous",
                                             `2` = "Homozygous Kaga",
                                             "UNCLASSIFIED" = "Unclassified")) %>% 
  # plot
  ggplot(aes(STATE_RECODE, FREQ, colour = STATE, fill = STATE)) +
  geom_violin() +
  geom_boxplot(width = .3) +
  ggbeeswarm::geom_quasirandom(color="#7D8491", size=0.4, alpha=0.9) +
  theme_bw() +
  scale_colour_manual(values = pal_hom_het_2_lines) +
  scale_fill_manual(values = pal_hom_het_2) +
  guides(colour = "none", fill = "none") +
  xlab("Genotype") +
  ylab("Proportion of reference bases covered") +
  ggtitle(paste("Site filter: ", 
                SITE_FILTER %>% 
                  stringr::str_replace_all("_", " "),
                "\nBin length: ",
                BIN_LENGTH,
                sep = ""))

## Save
ggsave(BASE_COV_TOT,
       base_cov_tot,
       device = "png",
       width = 10,
       height = 5.8,
       units = "in",
       dpi = 400)

# By chromosome

## Process
base_cov_df_chr = df %>% 
  split(., f = .$LANE) %>% 
  purrr::map(., function(LANE){
    # convert to ranges object
    lane_ranges = GenomicRanges::makeGRangesFromDataFrame(LANE,
                                                          keep.extra.columns = T,
                                                          ignore.strand = T,
                                                          seqnames.field = "chr", 
                                                          start.field = "BIN_START",
                                                          end.field = "BIN_END")
    # loop over each chromosome
    purrr::map(med_chr_lens$chr, function(CHR){
      # get total length of target chromosome
      target_chr_len = med_chr_lens$end[med_chr_lens$chr == CHR]
      # get total bases covered by each state per chromosome
      purrr::map_dfr(states, function(STATE){
        lane_ranges[lane_ranges$state == STATE & lane_ranges@seqnames == CHR] %>% 
          # merge contiguous ranges
          GenomicRanges::reduce(.) %>% 
          # get width of ranges
          width(.) %>% 
          # get total bases covered
          sum(.) %>% 
          # coerce into data frame
          data.frame("BASES_COVERED" = .)
      }, .id = "STATE") %>% 
        # add FREQ column
        dplyr::mutate(FREQ = BASES_COVERED / target_chr_len ) %>% 
        # add UNCLASSIFIED row
        tibble::add_row(STATE = "UNCLASSIFIED", 
                        BASES_COVERED = target_chr_len - sum(.$BASES_COVERED),
                        FREQ = (target_chr_len - sum(.$BASES_COVERED)) / target_chr_len)
    }) %>% 
      dplyr::bind_rows(.id = "CHR")
  }
  ) %>% 
  dplyr::bind_rows(.id = "LANE")

## Plot
base_cov_chr = base_cov_df_chr %>% 
  dplyr::mutate(STATE = factor(STATE, levels = c(0,1,2, "UNCLASSIFIED")),
                STATE_RECODE = dplyr::recode(STATE,
                                             `0` = "Homozygous Cab",
                                             `1` = "Heterozygous",
                                             `2` = "Homozygous Kaga",
                                             "UNCLASSIFIED" = "Unclassified"),
                CHR = factor(CHR, levels = med_chr_lens$chr)) %>% 
  # plot
  ggplot(aes(STATE_RECODE, FREQ, colour = STATE, fill = STATE)) +
  geom_violin() +
  geom_boxplot(width = .1) +
  ggbeeswarm::geom_quasirandom(color="#7D8491", size=0.1, alpha=0.7) +
  theme_bw() +
  scale_colour_manual(values = pal_hom_het_2_lines) +
  scale_fill_manual(values = pal_hom_het_2) +
  guides(colour = "none", fill = "none") +
  xlab("Genotype") +
  ylab("Proportion of reference bases covered") +
  facet_wrap(~CHR, nrow = 4, ncol = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste("Site filter: ", 
                SITE_FILTER %>% 
                  stringr::str_replace_all("_", " "),
                "\nBin length: ",
                BIN_LENGTH,
                sep = ""))

## Save
ggsave(BASE_COV_BY_CHROM,
       base_cov_chr,
       device = "png",
       width = 16,
       height = 13,
       units = "in",
       dpi = 400)

#######################
# Total sites covered by each state
#######################

# All sites

prop_sites_tot = df %>% 
  # get counts of sites per LANE and state
  dplyr::group_by(LANE, state) %>% 
  dplyr::count() %>% 
  # spread to one row per LANE
  tidyr::pivot_wider(id_cols = LANE, names_from = state, values_from = n) %>% 
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
  dplyr::mutate(STATE = factor(STATE, levels = c(0,1,2)),
                STATE_RECODE = dplyr::recode(STATE,
                                             `0` = "Homozygous Cab",
                                             `1` = "Heterozygous",
                                             `2` = "Homozygous Kaga")) %>% 
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
  ggtitle(paste("Site filter: ", 
                SITE_FILTER %>% 
                  stringr::str_replace_all("_", " "),
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

# Per chromosome

prop_sites_by_chrom = df %>% 
  dplyr::mutate(state = factor(state, levels = 0:2)) %>% 
  # get counts of sites per LANE and state
  dplyr::group_by(LANE, chr, state) %>%
  dplyr::count(.drop = F) %>% 
  # spread to one row per LANE
  tidyr::pivot_wider(id_cols = c(LANE, chr), names_from = state, values_from = n) %>% 
  # replace NAs with 0 manually , because `.drop = F` in `count` above doesn't work 
  dplyr::mutate(dplyr::across(c(`0`, `1`, `2`),
                              ~tidyr::replace_na(.x, 0))) %>% 
  # calculate frequencies of states per LANE
  dplyr::mutate(TOTAL = sum(`0`, `1`, `2`, na.rm = T),
                FREQ_0 = `0` / TOTAL,
                FREQ_1 = `1` / TOTAL,
                FREQ_2 = `2` / TOTAL) %>% 
  # gather
  tidyr::pivot_longer(cols = starts_with("FREQ_"),
                      names_to = "STATE",
                      names_prefix = "FREQ_",
                      values_to = "FREQ") %>% 
  # order STATE and recode with meaning
  dplyr::mutate(STATE = factor(STATE, levels = c(0,1,2)),
                STATE_RECODE = dplyr::recode(STATE,
                                             `0` = "Homozygous Cab",
                                             `1` = "Heterozygous",
                                             `2` = "Homozygous Kaga")) %>% 
  # plot
  ggplot(aes(STATE_RECODE, FREQ, colour = STATE, fill = STATE)) +
  geom_violin() +
  geom_boxplot(width = .1) +
  ggbeeswarm::geom_quasirandom(color="#7D8491", size=0.1, alpha=0.7) +
  theme_bw() +
  scale_colour_manual(values = pal_hom_het_2_lines) +
  scale_fill_manual(values = pal_hom_het_2) +
  guides(colour = "none", fill = "none") +
  xlab("Genotype") +
  ylab("Frequency") +
  facet_wrap(~chr, nrow = 4, ncol = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste("Site filter: ", 
                SITE_FILTER %>% 
                  stringr::str_replace_all("_", " "),
                "\nBin length: ",
                BIN_LENGTH,
                sep = ""))

ggsave(PROP_SITES_BY_CHROM,
       prop_sites_by_chrom,
       device = "png",
       width = 16,
       height = 13,
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
  split(., f = .$LANE) %>% 
  purrr::map(., function(LANE){
    # loop over CHR
    LANE %>% 
      split(., f = .$chr) %>% 
      purrr::map(., function(CHR){
        # Get lengths of each contiguous state
        cont_len = rle(CHR$state)
        
        # Get cumulative sum of those lengths
        cum_blocks = cumsum(cont_len$lengths)
        
        # Get rows that correspond to block changes
        block_bounds = CHR[cum_blocks, ] %>% 
          # Add end of previous block
          dplyr::mutate(END_PREV = dplyr::lag(BIN_END)) %>% 
          # Replace the NA in the first row with `1`
          dplyr::mutate(END_PREV = tidyr::replace_na(END_PREV, 1)) %>% 
          # Add colour
          dplyr::mutate(COLOUR = dplyr::recode(state,
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
                            paste("Site filter: ", 
                                  SITE_FILTER %>% 
                                    stringr::str_replace_all("_", " "),
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=4)

# Add rectangles in loop
counter = 0
purrr::map(block_bounds_list, function(LANE){
  # Add to counter_B
  counter <<- counter + 1
  # Add rectangles
  karyoploteR::kpRect(kp,
                      chr = LANE$chr,
                      x0 = LANE$END_PREV,
                      x1 = LANE$BIN_END,
                      y0 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(lower),
                      y1 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(upper),
                      col = LANE$COLOUR,
                      border = NA)
  # Add axis label
  karyoploteR::kpAddLabels(kp, labels = unique(LANE$LANE),
                           r0 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(lower),
                           r1 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(upper),
                           cex = 0.5)
})


dev.off()  


# Plot karyoplot WITH missing blocks

## Process
block_bounds_list = df %>% 
  # loop over LANE
  split(., f = .$LANE) %>% 
  purrr::map(., function(LANE){
  
    STRAIN = unique(LANE$LANE)
    # Create list of possible bins
    poss_bins = purrr::map(med_chr_lens$chr, function(CHR){
      # Get chr end
      CHR_END = med_chr_lens %>% 
        dplyr::filter(chr == CHR) %>% 
        dplyr::pull(end) %>% 
        as.numeric()
      # Get bin starts
      out = tibble::tibble(chr = as.numeric(CHR),
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
                              LANE %>% 
                                dplyr::select(chr, BIN_START, BIN_END, state),
                              by = c("chr", "BIN_START", "BIN_END")) %>% 
      # replace NAs with `UNCLASSIFIED`
      dplyr::mutate(state = state %>% 
                      tidyr::replace_na("UNCLASSIFIED"),
                    # add STRAIN
                    LANE = STRAIN) %>% 
      # add COLOUR
      dplyr::mutate(COLOUR = dplyr::recode(state,
                                           !!!pal_hom_het_2))
  
            
  })

# Plot karyoplot

png(file=KARYOPLOT_WIMISS,
    width=7800,
    height=23400,
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Site filter: ", 
                                  SITE_FILTER %>% 
                                    stringr::str_replace_all("_", " "),
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=4)

# Add rectangles in loop
counter = 0
purrr::map(block_bounds_list, function(LANE){
  # Add to counter
  counter <<- counter + 1
  # Add rectangles
  karyoploteR::kpRect(kp,
                      chr = LANE$chr,
                      x0 = LANE$BIN_START,
                      x1 = LANE$BIN_END,
                      y0 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(lower),
                      y1 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(upper),
                      col = LANE$COLOUR,
                      border = NA)
  # Add axis label
  karyoploteR::kpAddLabels(kp, labels = unique(LANE$LANE),
                           r0 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(lower),
                           r1 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(upper),
                           cex = 0.5)
})


dev.off()  
