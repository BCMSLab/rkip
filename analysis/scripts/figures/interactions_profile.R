# load libraries
library(tidyverse)
library(reshape2)
library(WGCNA)

# load data
load('data/rkip_wgcna.rda')

# global variables
figures_dir <- 'manuscript/figures'

# generate figure

ind <- c('PEBP1', 'PIK3C3', 'PIK3CB', 'TBC1D25', 'TBC1D5', 'TOLLIP', 'WDR45', 'WIPI1', 'TGFBR1')

expr <- mat[ind,] %>%
  melt %>%
  setNames(c('gene', 'sample', 'expr')) %>%
  left_join(design) %>%
  mutate(disease = factor(disease, levels = unique(disease))) %>%
  group_by(gene, disease) %>%
  mutate(ave = mean(expr), sd = sd(expr), low = ave - sd, upper = ave + sd)

corr <- cor(t(mat[ind,])) %>%
  melt %>%
  filter(Var1 == 'PEBP1') %>%
  dplyr::select(Var2, value) %>%
  setNames(c('gene', 'corr'))

pval <- corPvalueStudent(cor(t(mat[ind,])), ncol(mat)) %>%
  melt %>%
  filter(Var1 == 'PEBP1') %>%
  dplyr::select(Var2, value) %>%
  setNames(c('gene', 'pval'))

weights <- exportNetworkToVisANT(net$adj, threshold = .1) %>%
  filter(from == 'PEBP1') %>%
  dplyr::select(to, weight) %>%
  setNames(c('gene', 'weight'))

(left_join(expr, corr) %>%
  left_join(pval) %>%
  left_join(weights) %>%
  filter(gene != 'PEBP1') %>%
  mutate(pval = ifelse(pval < .01, '< 0.01', as.character(round(pval, 2)))) %>%
  ggplot(aes(x = '', y = expr, color = disease, group = disease)) +
  geom_point(position=position_dodge(width=0.7), aes(y=ave)) +
  geom_errorbar(position=position_dodge(width=0.7), aes(ymin = low, ymax = upper)) +
  geom_text(aes(y = 3, label = round(corr, 2)), color = 'black') +
  geom_text(aes(y = 2, label = pval), color = 'black') +
  scale_y_continuous(limits = c(2,14),
                     breaks = c(2, 3,5,7.5,10,12.5),
                     labels = c('p-value', "Pearson's", 5,7.5,10,12.5),
                     name = 'Log expression') +
  facet_wrap(~gene, nrow = 1) +
  theme_bw() +
  theme(legend.position = 'top') +
  labs(x = '', color = '')) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'interactions_profile.png', sep = '/'),
         width = 20, height = 10, units = 'cm')
