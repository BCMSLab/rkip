# load required libraries
library(tidyverse)
library(WGCNA)
library(xtable)

# load data
load('data/rkip_wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table
df <- exportNetworkToVisANT(net$adj, threshold = .1)
string <- read_tsv('data/string_interactions.tsv')

full_join(df, string) %>%
  filter(from == 'PEBP1') %>%
  left_join(select(ann, symbol, category), by = c('to'='symbol')) %>%
  mutate(evidence = ifelse(is.na(evidence), 'novel', evidence)) %>%
  filter(evidence != 'combined_score') %>%
  mutate(evidence = str_split(evidence, '_', simplify = TRUE)[, 1]) %>%
  group_by(category, evidence) %>%
  summarise(gene = paste(unique(to), collapse = ', ')) %>%
  ungroup() %>%
  mutate(category = ifelse(duplicated(category), '', category)) %>%
  setNames(c('Category', 'Evidence', 'Gene')) %>%
  xtable(caption = 'Summary of previously reported and novel PEBP1 interactions.',
         align = 'clll',
         label = 'tab:previousley_reported') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        sanitize.text.function = identity,
        comment = FALSE,
        file = paste(tables_dir, 'previously_reported.tex', sep = '/'))
