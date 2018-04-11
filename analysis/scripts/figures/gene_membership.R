# load libraries
library(tidyverse)
library(reshape2)
library(WGCNA)

# load data
load('data/rkip_wgcna.rda')

# golobal variables
figures_dir = 'manuscript/figures'

# generate figure
png(paste(figures_dir, 'gene_membership.png', sep = '/'),
    width = 10, height = 10, units = 'cm',
    res = 600, pointsize = 1)

net$diss %>%
  melt %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  acast(Var1 ~ Var2) %>%
  TOMplot(dendro = hclust(as.dist(net$diss), method = 'average'),
          Colors = net$merged_colors)

dev.off()
