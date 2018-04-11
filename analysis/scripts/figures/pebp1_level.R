# load libraries
library(tidyverse)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir <- 'manuscript/figures'

# generate figure
df <- read_tsv('data/atlas_expression.tsv')

(df %>%
  filter(`Gene Name` %in% c('GAPDH', 'PEBP1')) %>%
  select(`Gene Name`, starts_with('GSM')) %>%
  group_by(`Gene Name`) %>%
  summarise_all(mean) %>%
  gather(sample, expr, -`Gene Name`) %>%
  left_join(design) %>%
  mutate(disease = factor(disease, levels = unique(disease))) %>%
  spread(`Gene Name`, expr) %>%
  mutate(fc = PEBP1 - GAPDH) %>%
  group_by(disease) %>%
  summarise(fc  = mean(fc)) %>%
  mutate(fc = fc/filter(., disease == 'benign prostate tumor') %>% pull(fc)) %>%
  ggplot(aes(x = disease, y = fc)) +
  geom_col(width = .7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Relative fold-change\n')) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'pebp1_level.png', sep = '/'),
         width = 7, height = 10, units = 'cm')
