# load libraries
library(tidyverse)
library(reshape2)
library(ggimage)
library(cowplot)

# load data
if(!file.exists('data/image.png')) {
  download.file('https://ndownloader.figshare.com/files/12293744',
                destfile = 'data/image.png')
}

g1 <- data_frame(image = 'data/image.png') %>%
  ggplot() +
  geom_image(aes(x = .5, y = .5, image = image), size = 1) +
  scale_x_continuous(breaks = c(0.125, 0.375, 0.625, 0.875),
                     labels = c('', 'RKIP', 'Hoechst', 'Merged'),
                     position = 'top',
                     name = NULL) +
  scale_y_continuous(breaks = c(.083, .25, 0.416, 0.583, .75 , 0.916),
                     labels = c('MAP1LC3B', 'WIPI1', 'TBC1D5', 'TOLLIP', 'PIK3CB', 'PIK3C3'),
                     name = NULL) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 8),
        axis.line = element_blank())

if(!file.exists('data/coloc.csv')) {
  download.file('https://ndownloader.figshare.com/files/12293747',
                destfile = 'data/coloc.csv')
}
df <- read_csv('data/coloc.csv') %>%
  setNames(c('symbol', "Pearson's", "Manders (M1)", "Manders (M2)"))

g2 <- df %>%
  gather(type, value, -symbol) %>%
  ggplot(aes(x = symbol, y = value)) +
  geom_jitter(width = .3, alpha = .5) +
  geom_point(data = df %>%
                  gather(type, value, -symbol) %>%
                  group_by(symbol, type) %>%
                  summarise(ave = mean(value)),
                aes(y = ave),
                color = 'red',
                width = .3) +
  geom_errorbar(data = df %>%
                gather(type, value, -symbol) %>%
                group_by(symbol, type) %>%
                summarise(ave = mean(value),
                          sd = sd(value),
                          upper = ave + sd,
                          lower = ave - sd),
              aes(y = ave, ymin = lower, ymax = upper),
              color = 'red',
              width = .3) +
  facet_wrap(~type, nrow = 1) +
  lims(y = c(0,1.05)) +
  labs(y = 'Coefficient Value\n', x = '') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_grid(g1, g2,
          ncol = 1, rel_heights = c(2,1),
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10,
          scale = .95) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/colocalization.png',
         width = 18, height = 25, units = 'cm')
