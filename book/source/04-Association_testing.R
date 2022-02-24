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

run_gwas <- function(d,m,p,invers_norm=F) {
  ids = paste("snp", 1:nrow(m), sep="")
  mm = data.frame(ids, m)
  colnames(d) = ids
  colnames(mm) = c("snp", "Chr", "pos")
  mm = mm[,1:3]
  map = mm
  n_geno=nrow(d)
  Plots = matrix(1:1200,nc = 10)
  Plot_Row = row(Plots)
  Plot_Col = col(Plots)
  data = data.frame(Geno = paste0('Geno',1: n_geno), Plot = sample(Plots)[1:nrow(d)])
  data$Row = Plot_Row[data$Plot]
  data$Col = Plot_Col[data$Plot]
  if(invers_norm) {
    data$y = my.invnorm(p[,2])
  } else {
    data$y = p[,2]
  }
  X = as.matrix(d)
  X[is.na(X)]=0
  row.names(X) = data[,1]
  X_centered = sweep(X,2,colMeans(X),'-') # center marker genotypes
  K = tcrossprod(X_centered) / ncol(X_centered)
  rownames(K) = colnames(K) = data$Plot
  field = data[,c('Row','Col')]
  dists = as.matrix(dist(field))
  h = median(dists)
  K_plot = gausskernel(field,h^2/2); diag(K_plot)=1 # 
  rownames(K_plot) = colnames(K_plot) = data$Plot
  gwas = GridLMM_GWAS(
    formula = y~1 + (1|Geno) + (1|Plot), # the same error model is used for each marker. It is specified similarly to lmer
    test_formula = ~1, # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
    reduced_formula = ~0, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
    data = data, # The dataframe to look for terms from the 3 models
    weights = NULL, # optional observation-specific weights
    X = X, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
    X_ID = 'Geno', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
    h2_start = NULL, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using GridLMM_ML
    h2_step = 0.01, # step size per random effect for testing alternate values of h2
    max_steps = 100, # maximum number of steps of size h2_step to take from h2_start
    X_map = map, # Optional. The marker positions.
    relmat = list(Plot = K), # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
    centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
    scaleX = FALSE, # Should the markers be scaled to have constant variance when calculating the GRM?
    fillNAX = FALSE, # Should missing marker data be filled in with the mean allele frequency?
    method = 'REML', # REML = Wald test, ML = LRT, BF = calculate Bayes factors
    mc.cores = my_detectCores(), # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
    verbose = FALSE # Should progress be printed to the screen?
  )
  return(gwas)
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

