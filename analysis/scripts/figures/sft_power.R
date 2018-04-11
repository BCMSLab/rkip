# load libraries
library(tidyverse)
library(reshape2)
library(cowplot)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir = 'manuscript/figures'


# generate figure

(plot_grid(sft$fitIndices %>%
             mutate(y = -sign(slope)*SFT.R.sq) %>%
             ggplot(aes(x = Power, y = y)) +
             geom_point() +
             geom_hline(yintercept = 0.39, lty = 2, color = 'red') +
             labs(x = 'Soft Thresold (power)',
                  y = 'Scale-Free Topology') +
             theme_bw(),
           sft$fitIndices %>%
             ggplot(aes(x = Power, y = `mean.k.`)) +
             geom_point() +
             geom_vline(xintercept = 5, lty = 2, color = 'red') +
             labs(x = 'Soft Threshold (power)',
                  y = 'Mean Connectivity') +
             theme_bw(),
           labels = 'AUTO',
           label_size = 10,
           label_fontface = 'plain',
           nrow = 1,
           scale = .9)) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'sft_power.png', sep = '/'),
         width = 16, height = 7, units = 'cm')
