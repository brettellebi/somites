########################
# Libraries
########################

library(tidyverse)
library(PhenotypeSimulator)
library(GridLMM)
library(DT)
library(KRLS)
library(cowplot)
library(plotly)
library(gtools)
library(ggbeeswarm)

########################
# Paths
########################

working_dir = "/hps/nobackup/birney/users/ian/somites"
lts_dir = "/nfs/research/birney/users/ian/somites"

########################
# GridLMM functions
########################

my.invnorm = function(x) {
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}

########################
# Generic functions
########################

clean_gwas_res = function(gwas_results, bin_length, chr_lens){
  gwas_results$results %>% 
    dplyr::left_join(med_chr_lens, by = c("Chr" = "chr")) %>% 
    # add x-coord
    dplyr::mutate(X_COORD = pos + TOT) %>% 
    # change column names
    dplyr::rename(CHROM = Chr, BIN_START = pos) %>% 
    # add BIN_END
    dplyr::mutate(BIN_END = BIN_START + bin_length - 1) %>% 
    # add locus
    dplyr::mutate(LOCUS = paste(CHROM, BIN_START, sep = ":")) 
}

plot_man = function(df, site_filter, phenotype, bin_length, gwas_pal, size = 0.5, alpha = 0.5, med_chr_lens, sig_level = NULL){
  # Create palette
  pal = rep_len(gwas_pal, length.out = nrow(med_chr_lens))
  names(pal) = med_chr_lens$chr
  
  df = df %>% 
    # create `COLOUR` vector
    dplyr::mutate(COLOUR = dplyr::case_when(!is.null(sig_level) & p_value_REML < sig_level ~ gwas_pal[1],
                                            gtools::even(CHROM) ~ gwas_pal[2],
                                            gtools::odd(CHROM) ~ gwas_pal[3])) %>% 
    dplyr::mutate(CHROM = factor(CHROM, levels = med_chr_lens$chr)) 
  
  out_plot = df %>% 
    ggplot(aes(x = X_COORD,
               y = -log10(p_value_REML),
               label = BIN_START,
               label2 = BIN_END)) + 
    geom_point(colour = df$COLOUR,
               size = size,
               alpha = alpha) +
    #scale_color_manual(values = gwas_pal) +
    scale_x_continuous(breaks = med_chr_lens$MID_TOT, 
                       labels = med_chr_lens$chr) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) +
    guides(colour = "none") +
    ggtitle(paste("Site filter: ", site_filter, "\nPhenotype: ", phenotype, "\nBin length: ",  bin_length, sep = "")) +
    xlab("Chromosome") +
    ylab("-log10(p-value)") + 
    geom_hline(yintercept = -log10(sig_level), colour = "#60D394", linetype = "dashed")
  
  return(out_plot)
  
}

plot_int_man = function(df, phenotype, bin_length, gwas_pal, size = 0.5, alpha = 0.5, med_chr_lens, sig_line = NULL){
  # Create palette
  pal = rep_len(gwas_pal, length.out = nrow(med_chr_lens))
  names(pal) = med_chr_lens$chr
  
  # Create plot
  p = df %>% 
    dplyr::mutate(CHROM = factor(CHROM, levels = med_chr_lens$chr)) %>% 
    ggplot(aes(x = X_COORD,
               y = -log10(p_value_REML),
               colour = CHROM,
               label = BIN_START,
               label2 = BIN_END)) + 
    geom_point(size = size,
               alpha = alpha) +
    scale_color_manual(values = pal) +
    scale_x_continuous(breaks = med_chr_lens$MID_TOT, 
                       labels = med_chr_lens$chr) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) +
    guides(colour = "none") +
    ggtitle(paste("Phenotype: ", phenotype, "\nBin length: ",  bin_length, sep = "")) +
    xlab("Chromosome") +
    ylab("-log10(p-value)") + 
    geom_hline(yintercept = sig_line, colour = "#1effbc", linetype = "dashed") 
  
  ggplotly(p, tooltip = c("CHROM", "BIN_START", "BIN_END"))
}

########################
# Plotting parameters
########################

gwas_pal = c("#2B2D42", "#F7B267", "#F25C54")
names(gwas_pal) = c("target", "even chr", "odd chr")
significance_line = 3.6
suggestive_line = 2.9

# Intercept
intercept_pal = c("#EF476F", "#8D99AE", "#2b2d42")
names(intercept_pal) = c("target", "even chr", "odd chr")

# Mean
mean_pal = c("#D81E5B", "#8AA399", "#084C61")
names(mean_pal) = c("target", "even chr", "odd chr")

# PSM
unsegmented_psm_area_pal = c("#E59500", "#D9D0DE", "#401F3E")
names(mean_pal) = c("target", "even chr", "odd chr")

########################
# HdrR chromosome data
########################
# Get chromosome lengths
med_chr_lens = read.table(here::here("data",
                                     "Oryzias_latipes.ASM223467v1.dna.toplevel.fa_chr_counts.txt"),
                          col.names = c("chr", "end"))
# Add start
med_chr_lens$start = 1
# Reorder
med_chr_lens = med_chr_lens %>% 
  dplyr::select(chr, start, end) %>% 
  # remove MT
  dplyr::filter(chr != "MT") %>% 
  # convert to integer
  dplyr::mutate(chr = as.integer(chr)) %>% 
  # Add cumulative bases
  dplyr::mutate(CUMSUM = cumsum(end),
                TOT = CUMSUM - end) %>% 
  # Add midpoint for each chr
  dplyr::mutate(MID_TOT = TOT + (end / 2))

