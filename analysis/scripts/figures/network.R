# load libraries
library(tidyverse)
library(ggraph)
library(igraph)
library(RSQLite)
library(cRegulome)
library(WGCNA)
library(reshape2)
library(tidygraph)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir = 'manuscript/figures'

# generate figure

ind <- c('PEBP1', 'PIK3C3', 'PIK3CB', 'TBC1D25', 'TBC1D5', 'TOLLIP', 'WDR45', 'WIPI1', 'TGFBR1')

conn <- dbConnect(SQLite(), 'data/cRegulome.db')

tf <- get_tf(conn,
             tf = c('ERCC6', 'VEZF1'),
             study = 'PRAD',
             targets_only = TRUE) %>%
  filter(feature %in% ind) %>%
  dplyr::select(-study) %>%
  setNames(c('from', 'to', 'cor'))

mir <- get_mir(conn,
               mir = c('hsa-mir-23c', 'hsa-mir-378c', 'hsa-mir-761'),
               study = 'PRAD') %>%
  filter(feature %in% ind) %>%
  dplyr::select(-study) %>%
  setNames(c('from', 'to', 'cor'))

dbDisconnect(conn)

corr <- cor(t(mat[ind,]))

gi <- corr %>%
  melt %>%
  filter(abs(value) > .5) %>%
  filter(value != 1) %>%
  setNames(c('from', 'to', 'cor'))

node_type <- list(miRNA = unique(mir$from),
                   TF = unique(tf$from),
                   Gene = ind) %>%
  melt() %>%
  setNames(c('name', 'Type')) %>%
  mutate_all(as.character)

(list(GI = gi, TF = tf, miRNA = mir) %>%
  bind_rows(.id = 'Type') %>%
  dplyr::select(from, to, cor, Type) %>%
  mutate(Direction = ifelse(cor > 0, 'Positive', 'Negative')) %>%
  as_tbl_graph(dircted = FALSE) %>%
  left_join(node_type) %>%
  activate(edges) %>%
  mutate(inc = edge_is_incident(9)) %>%
  ggraph(layout = 'kk') +
  geom_edge_link(aes(color = Direction)) +
  geom_node_point(aes(color = Type), size = 5) +
  geom_node_text(aes(label = name)) +
  theme_graph() +
  theme(legend.direction = 'vertical')) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'network.png', sep = '/'),
         width = 15, height = 15, units = 'cm')
