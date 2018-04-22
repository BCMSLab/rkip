# load libraries
library(tidyverse)
library(cowplot)
library(RSQLite)
library(cRegulome)

# global vars
figures_dir = 'manuscript/figures'

# generate figure

ind <- c('PEBP1', 'PIK3C3', 'PIK3CB', 'TBC1D25', 'TBC1D5', 'TOLLIP', 'WDR45', 'WIPI1', 'TGFBR1')

conn <- dbConnect(SQLite(), 'data/cRegulome.db')

tf <- get_tf(conn,
             tf = c('ERCC6', 'VEZF1'),
             study = 'PRAD',
             targets_only = TRUE) %>%
  filter(feature %in% ind)

mir <- get_mir(conn,
               mir = c('hsa-mir-23c', 'hsa-mir-378c', 'hsa-mir-761'),
               study = 'PRAD') %>%
  filter(feature %in% ind)

dbDisconnect(conn)

(list(TF = tf,
      microRNA = mir) %>%
    map(function(x) {
      names(x)[1] <- 'reg'

      x %>%
        ggplot(aes(x = feature, y = cor)) +
        facet_grid(~reg, scales = 'free_x', space = 'free_x') +
        geom_col() +
        geom_hline(yintercept = 0, color = 'gray') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = '', y = 'Expression Correlation') +
        lims(y = c(-1, 1))
    })) %>%
  plot_grid(plotlist = .,
            rel_widths = c(1,2),
            scale = .9,
            labels = 'AUTO',
            label_fontface = 'plain',
            label_size = 10) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'common_regulation.png', sep = '/'),
         width = 20, height = 10, units = 'cm')
