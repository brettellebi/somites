# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN_F01 = here::here("data/F0_F1_period.xlsx")
IN_F2 = here::here("config/phenos_with_reporter_genoandpheno.csv")

## True
IN_F01 = snakemake@input[["f01"]]
IN_F2 = snakemake@input[["f2"]]
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]

########################
# Plotting parameters
########################

# Intercept
intercept_pal = c("#8D99AE", "#2b2d42")
#names(intercept_pal) = c("even chr", "odd chr")

# Mean
mean_pal = c("#8AA399", "#084C61")
#names(mean_pal) = c("even chr", "odd chr")

# PSM
unsegmented_psm_area_pal = c("#D9D0DE", "#401F3E")

# Get lighter/darker functions

devtools::source_gist("c5015ee666cdf8d9f7e25fa3c8063c99")

########################
# Read in file
########################

df_f2 = readr::read_delim(IN_F2, delim = ";") %>% 
  # add `GEN` column
  dplyr::mutate(GEN = "F2")


# Read in F0 and F1 data
df_f01 = readxl::read_xlsx(IN_F01) %>% 
  dplyr::mutate(sample = fish) %>% 
  dplyr::mutate(GEN = dplyr::case_when(str_detect(fish, "^C") ~ "F0",
                                       str_detect(fish, "^K") ~ "F1"))

# Bind two data frames
df_all = dplyr::bind_rows(df_f01, df_f2) %>% 
  # factorise Microscope
  dplyr::mutate(Microscope = factor(Microscope, levels = c("AU", "DB")))

########################
# Kruskal-Wallis test
########################

kw_df = df_all %>% 
  # pivot longer to put phenotypes values in one column
  tidyr::pivot_longer(cols = c(mean, intercept, unsegmented_psm_area),
                      names_to = "phenotype",
                      values_to = "value") %>% 
  dplyr::group_by(phenotype) %>% 
  tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data,
                                   ~kruskal.test(x = .$value, g = .$Microscope))) %>%
  dplyr::select(-data) %>% 
  dplyr::mutate(model_tidy = purrr::map(model, broom::tidy)) %>%
  tidyr::unnest(model_tidy) %>% 
  rstatix::add_significance(p.col = "p.value") %>% 
  # remove model
  dplyr::select(-model) %>% 
  # reduce to 3 digits
  dplyr::mutate(p.value = signif(p.value, digits = 3)) %>% 
  # paste p-value with significance
  dplyr::mutate(p_final = dplyr::case_when(p.value.signif == "ns" ~ paste("p =", p.value),
                                           TRUE ~ paste("p =", p.value, p.value.signif))) %>% 
  # add `Microscope` column with 'DB' so that the text maps there on the plots
  dplyr::mutate(Microscope = factor("DB", levels = c("AU", "DB")))

########################
# Plot
########################

########### Mean

mean_fig = df_all %>% 
  # remove NAs
  dplyr::filter(!is.na(Microscope)) %>% 
  ggplot(aes(GEN, mean, fill = Microscope)) +
  geom_violin() + 
  geom_boxplot(width = 0.3) +
  ggbeeswarm::geom_beeswarm(aes(GEN, mean, colour = Microscope), size = 0.4, alpha = 0.5) +
  facet_wrap(~Microscope, ncol = 2) + 
  scale_colour_manual(values = lighter(c("#8D99AE", "#2b2d42"), amount = 50)) +
  scale_fill_manual(values = darker(c("#8D99AE", "#2b2d42"), amount = 50)) +
  cowplot::theme_cowplot() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_text(face = "bold")) +
  xlab("generation") +
  ylab("mean period") +
  guides(fill = "none",
         colour = "none") + 
  # add p-value
  geom_text(data = kw_df %>% 
                       dplyr::filter(phenotype == "mean"),
                     aes(x = "F2", y = -Inf, label = p_final,
                         vjust = -1
                         ))

########### Intercept

intercept_fig = df_all %>% 
  # remove NAs
  dplyr::filter(!is.na(Microscope)) %>% 
  ggplot(aes(GEN, intercept, fill = Microscope)) +
  geom_violin() + 
  geom_boxplot(width = 0.3) +
  ggbeeswarm::geom_beeswarm(aes(GEN, intercept, colour = Microscope), size = 0.4, alpha = 0.5) +
  facet_grid(cols = vars(Microscope)) + 
  scale_colour_manual(values = lighter(c("#177e89", "#084c61"), amount = 50)) +
  scale_fill_manual(values = darker(c("#177e89", "#084c61"), amount = 50)) +
  cowplot::theme_cowplot() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_text(face = "bold")) +
  xlab("generation") +
  ylab("period intercept") +
  guides(fill = "none",
         colour = "none") +
  # add p-value
  geom_text(data = kw_df %>% 
              dplyr::filter(phenotype == "intercept"),
            aes(x = "F2", y = -Inf, label = p_final,
                vjust = -1
            ))

########### PSM

psm_fig = df_all %>% 
  # remove NAs
  dplyr::filter(!is.na(Microscope)) %>% 
  ggplot(aes(GEN, unsegmented_psm_area, fill = Microscope)) +
  geom_violin() + 
  geom_boxplot(width = 0.3) +
  ggbeeswarm::geom_beeswarm(aes(GEN, unsegmented_psm_area, colour = Microscope), size = 0.4, alpha = 0.5) +
  facet_grid(cols = vars(Microscope)) + 
  scale_colour_manual(values = lighter(c("#D9D0DE", "#401F3E"), amount = 50)) +
  scale_fill_manual(values = darker(c("#D9D0DE", "#401F3E"), amount = 50)) +
  cowplot::theme_cowplot() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_text(face = "bold")) +
  xlab("generation") +
  ylab("unsegmented PSM area") +
  guides(fill = "none",
         colour = "none") +
  # add p-value
  geom_text(data = kw_df %>% 
              dplyr::filter(phenotype == "unsegmented_psm_area"),
            aes(x = "F2", y = -Inf, label = p_final,
                vjust = -1
            ))


########### Together

period_final = cowplot::plot_grid(intercept_fig,
                                  mean_fig,
                                  psm_fig,
                                  align = "hv",
                                  nrow = 3,
                                  labels = c("A", "B", "C"),
                                  label_size = 16)

ggsave(OUT_PNG,
       device = "png",
       width = 11,
       height = 13.5,
       units = "in",
       dpi = 400)

ggsave(OUT_PDF,
       device = "pdf",
       width = 11,
       height = 13.5,
       units = "in",
       dpi = 400)


