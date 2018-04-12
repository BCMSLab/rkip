# load libraries
library(WGCNA)
library(xtable)
library(tidyverse)

# load data
load('data/rkip_wgcna.rda')

# global variables
tables_dir = 'manuscript/tables'

# generate table
exportNetworkToVisANT(net$adj, threshold = .1) %>%
  filter(from %in% split(ann$symbol, ann$category)$pebp) %>%
  left_join(genes, by = c('to'='symbol')) %>%
  left_join(ann, by = c('to'='symbol')) %>%
  filter(category != 'pebp') %>%
  group_by(color, from, category) %>%
  summarise(gene = paste(unique(to), collapse = ', ')) %>%
  spread(category, gene) %>%
  ungroup() %>%
  mutate(color = ifelse(duplicated(color), '', color)) %>%
  xtable(caption = 'PEB interactions with autophagy and EMT.',
         label = 'tab:peb_interactions',
         align = 'cllp{.4\\textwidth}p{.3\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        sanitize.text.function = identity,
        comment = FALSE,
        file = paste(tables_dir, 'peb_interactions.tex', sep = '/'))

