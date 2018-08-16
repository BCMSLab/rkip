# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

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
  left_join(genes) %>%
  na.omit() %>%
  mutate(symbol = paste('\\emph{', symbol, '}', sep = '')) %>%
  group_by(color, term) %>%
  summarise(gene = paste(unique(symbol), collapse = ', ')) %>%
  spread(term, gene) %>%
  setNames(c('Module', 'Autophagy',
             'Epithelial to Mesenchymal Transition',
             'Phosphatidylethanolamine Binding')) %>%
  xtable(caption = 'Gene members in different modules/colors.',
         label = 'tab:module_members',
         align = 'clp{.4\\textwidth}p{.2\\textwidth}p{.25\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        file = paste(tables_dir, 'module_members.tex', sep = '/'))
