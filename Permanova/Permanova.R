# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)
library(microbiomeMarker)

#import data

NBA <- import_qiime2(
  otu_qza = "table.qza", sam_tab = "preterm_metadata_bacteria_nutrition_updated_v3.txt")

detach("package:microbiomeMarker", unload = TRUE)  #unload microbiome maker again as it interferes with other commands later on

#rename
pseq <- NBA

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Nutrition_bifidobacterium", size = 3)
print(p)

# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ Nutrition_bifidobacterium,
                    data = meta, permutations=999, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["Nutrition_bifidobacterium", "Pr(>F)"])
# R2
print(as.data.frame(permanova$aov.tab)["Nutrition_bifidobacterium", "R2"])
