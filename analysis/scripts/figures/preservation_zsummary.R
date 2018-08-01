# load libraries
library(tidyverse)
library(reshape2)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir = 'manuscript/figures'

# generate figure

(module_preserv$quality$Z$ref.ref[-9] %>%
    map(rownames_to_column, var = 'color') %>%
    bind_rows(.id = 'study') %>%
    mutate(gse = str_split(study, '\\.', simplify = TRUE)[, 2],
           color = ifelse(color == 'gold', 'gray', color)) %>%
    ggplot(aes(x = moduleSize, y = Zsummary.qual, color = color)) +
    geom_point() +
    geom_abline(intercept = c(2,5), slope = 0, lty = 2) +
    theme_bw() +
    scale_color_manual(values = c('blue', 'brown', 'gray', 'yellow')) +
    theme(legend.position = 'top') +
    facet_wrap(~gse, nrow = 2) +
    labs(x = 'Module Size', y = 'Preservation Z Summary', color = '') +
    guides(color = guide_legend(nrow = 1))) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'preservation_zsummary.png', sep = '/'),
         width = 20, height = 12, units = 'cm')
