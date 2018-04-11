# load libraries
library(tidyverse)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir = 'manuscript/figures'

# generate figure
(cor(net$mes, as.numeric(as.factor(design$disease))) %>%
  as.data.frame() %>%
  mutate(color = c('blue', 'brown', 'yellow')) %>%
  ggplot(aes(x = color, y = V1)) +
  geom_col() +
  theme_bw() +
  lims(y = c(-.5, .5)) +
  labs(x = '', y = "Pearsons's Correlation") +
  geom_abline(intercept = 0, slope = 0, lty = 2)) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'module_cor.png', sep = '/'),
         width = 8, height = 7, units = 'cm')
