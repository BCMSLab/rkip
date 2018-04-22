# load libraries
library(tidyverse)
library(cowplot)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir = 'manuscript/figures'

# generate figure

plot_grid(
  cor(net$mes, as.numeric(as.factor(design$disease))) %>%
    as.data.frame() %>%
    mutate(color = c('blue', 'brown', 'yellow')) %>%
    ggplot(aes(x = color, y = V1)) +
    geom_col() +
    theme_bw() +
    lims(y = c(-.5, .5)) +
    labs(x = '', y = "Pearsons's Correlation") +
    geom_abline(intercept = 0, slope = 0, lty = 2),
  cmdscale(net$diss) %>%
    as.data.frame() %>%
    setNames(c('D1', 'D2')) %>%
    mutate(color = net$merged_colors) %>%
    ggplot(aes(x = D1, y = D2, color = color)) +
    geom_point() +
    scale_color_manual(values = c('blue', 'brown', 'yellow')) +
    theme_bw() +
    theme(legend.position = 'none'),
  labels = 'AUTO',
  label_size = 10,
  label_fontface = 'plain',
  scale = .9
) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'module_cor.png', sep = '/'),
         width = 16, height = 7, units = 'cm')
