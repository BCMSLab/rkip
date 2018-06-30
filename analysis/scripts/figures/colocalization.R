# load libraries
library(tidyverse)
library(reshape2)

# load data
download.file('https://ndownloader.figshare.com/files/12265169',
              destfile = 'data/coloc.csv')
df <- read_csv('data/coloc.csv')

# generate figure
pval <- with(df, split(pearsons, symbol)) %>%
  map(function(x) t.test(x)$p.value) %>%
  melt %>%
  setNames(c('pvalue', 'symbol'))

all(pval$pvalue < .001)

(df %>%
  ggplot(aes(x = symbol, y = pearsons)) +
  geom_jitter(width = .2, alpha = .7) +
  geom_point(data = group_by(df, symbol) %>%
               summarise(ave = mean(pearsons)),
             aes(y = ave),
             color = 'red') +
  geom_errorbar(data = group_by(df, symbol) %>%
               summarise(ave = mean(pearsons),
                         sd = sd(pearsons),
                         upper = ave + sd,
                         lower = ave - sd),
             aes(y = ave, ymin = lower, ymax = upper),
             color = 'red',
             width = .2) +
  geom_point(data = pval, aes(x = symbol, y = 1), pch = 8) +
  annotate('text', x = 10, y = 0, label = 'p-value < 0.001 *', hjust = 1) +
  labs(y = "Co-localization (Pearson's Coeff)", x = '') +
  lims(y = c(0, 1)) +
  theme_bw() +
  theme()) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/colocalization.png',
         width = 18, height = 8, units = 'cm')
