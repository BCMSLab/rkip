# load libraries
library(tidyverse)
library(xtable)

# load data
load('data/rkip_wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables'

# generate table
data_frame(term = c('epithelial to mesenchymal transition',
                    'autophagy',
                    'phosphatidylethanolamine binding'),
           go = c("GO:0001837", "GO:0006914", "GO:0008429")) %>%
  left_join(ann) %>%
  mutate(term = paste(term, ' (', go, ')', sep = '')) %>%
  group_by(term) %>%
  summarise(genes = paste(unique(symbol), collapse = ', ')) %>%
  setNames(c('Term', 'Genes')) %>%
  xtable(caption = 'Gene members of the three gene ontology terms.',
         align = 'clp{.6\\textwidth}',
         label = 'tab:go_genes') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        file = paste(tables_dir, 'go_genes.tex', sep = '/'))
