# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/rkip_wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables'

# generate table
map(md, function(x) data.frame(nsample = nrow(x$data),
                                  ngenes = ncol(x$data))) %>%
  bind_rows(.id = 'study') %>%
  filter(study != 'ref') %>%
  mutate(reference = c('Baca et al. Cell 2013',
                       'Barbieri et al. Nat Genet 2012',
                       'Kumar et al. Nat Med 2016',
                       'Gao et al., Cell 2014',
                       'Taylor et al. Cancer Cell 2010',
                       'Robinson et al. Cell. 2015',
                       'TCGA, Cell 2015',
                       'TCGA, Cell 2015')) %>%
  setNames(c('Study ID', 'Samples', 'Genes', 'Reference')) %>%
  xtable(caption = 'Studies of human prostate cancer subjects.',
         align = 'clccl',
         label = 'tab:datasets') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'datasets.tex', sep = '/'))
