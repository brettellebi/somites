#!/usr/bin/env Rscript

####################
# Libraries
####################

library(tidyverse)
library(karyoploteR)
library(ggbeeswarm)

####################
# Palettes
####################

pal_hom_het_1 = c("#947EB0", "#A3A5C3", "#A9D2D5", "#E1E1DF")
names(pal_hom_het_1) = c(0:2, "UNCLASSIFIED")

pal_hom_het_2 = c("#e28413", "#000022", "#DE3C4B", "#FBF5F3")
names(pal_hom_het_2) = c(0:2, "UNCLASSIFIED")

pal_hom_het_2_lines = c(karyoploteR::darker(pal_hom_het_2[1], 100),
                        karyoploteR::lighter(pal_hom_het_2[2], 100),
                        karyoploteR::darker(pal_hom_het_2[3], 100),
                        karyoploteR::darker(pal_hom_het_2[4]))
names(pal_hom_het_2_lines) = c(0:2, "UNCLASSIFIED")
