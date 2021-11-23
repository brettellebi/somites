#!/usr/bin/env Rscript

####################
# Libraries
####################

library(tidyverse)
library(karyoploteR)
library(ggbeeswarm)
#library(sommer)
library(plotly)
library(DT)

####################
# Palettes
####################

pal_hom_het_1 = c("#947EB0", "#A3A5C3", "#A9D2D5", "#E1E1DF")
names(pal_hom_het_1) = c(0:2, "UNCLASSIFIED")

pal_hom_het_2 = c("#43AA8B", "#000022", "#DE3C4B", "#FBF5F3")
names(pal_hom_het_2) = c(0:2, "UNCLASSIFIED")

pal_hom_het_2_lines = c(karyoploteR::darker(pal_hom_het_2[1], 100),
                        karyoploteR::lighter(pal_hom_het_2[2], 100),
                        karyoploteR::darker(pal_hom_het_2[3], 100),
                        karyoploteR::darker(pal_hom_het_2[4]))
names(pal_hom_het_2_lines) = c(0:2, "UNCLASSIFIED")

pal_ck_1 = c("#FF6F59", "#000022", "#43AA8B", "#FBF5F3")
names(pal_ck_1) = c(0:2, "UNCLASSIFIED")

pal_ck_1_lines = c(karyoploteR::darker(pal_ck_1[1], 100),
                   karyoploteR::lighter(pal_ck_1[2], 100),
                   karyoploteR::darker(pal_ck_1[3], 100),
                   karyoploteR::darker(pal_ck_1[4]))
names(pal_ck_1_lines) = c(0:2, "UNCLASSIFIED")
