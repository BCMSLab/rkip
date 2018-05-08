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
  mutate(study = gsub('_', '.', study)) %>%
  filter(study != 'ref') %>%
  mutate(reference = c('\\cite{Baca2013PunctuatedGenomes}',
                       '\\cite{Barbieri2012ExomeCancer}',
                       '\\cite{Kumar2016SubstantialCancer}',
                       '\\cite{Gao2014DeterministicNeocortex}',
                       '\\cite{Taylor2010IntegrativeCancer}',
                       '\\cite{Robinson2015IntegrativeCancer}',
                       '\\cite{CancerGenomeAtlasResearchNetwork2015TheCancer.}',
                       '\\cite{CancerGenomeAtlasResearchNetwork2015TheCancer.}')) %>%
  setNames(c('Study ID', 'Samples', 'Genes', 'Reference')) %>%
  xtable(caption = 'Studies of human prostate cancer subjects.',
         align = 'clccc',
         label = 'tab:datasets') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        sanitize.text.function = identity,
        comment = FALSE,
        file = paste(tables_dir, 'datasets.tex', sep = '/'))
