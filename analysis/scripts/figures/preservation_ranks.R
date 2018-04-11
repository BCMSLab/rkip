# load libraries
library(tidyverse)
library(reshape2)
library(WGCNA)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir = 'manuscript/figures'

# generate figure

(module_preserv$quality$observed$ref.ref[-9] %>%
    map(rownames_to_column, var = 'color') %>%
    bind_rows(.id = 'study') %>%
    mutate(gse = str_split(study, '\\.', simplify = TRUE)[, 2]) %>%
    ggplot(aes(x = moduleSize, y = medianRank.qual, color = color)) +
    geom_point() +
    facet_wrap(~study, nrow = 2) +
    theme_bw() +
    scale_color_manual(values = c('blue', 'brown', 'gold', 'yellow')) +
    theme(legend.position = 'top') +
    labs(x = 'Module Size', y = 'Preservation Median Rank', color = '') +
    guides(color = guide_legend(nrow = 1))) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'preservation_ranks.png', sep = '/'),
         width = 20, height = 12, units = 'cm')
