# load required libraries
library(tidyverse)
library(reshape2)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(trackViewer)
library(cowplot)

# global variables
figures_dir <- 'manuscript/figures'

# generate figure
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

gene_ids <- list(PEBP1 = 5037, TBC1D5 = 9779, FLYWCH2 = 114984,
                 PIK3C3 = 5289, WDR45 = 11152, TMEM230 = 29058)

gene_coords <- map(gene_ids, function(x) {
  gr <- genes(txdb, filter = list(gene_id = x))
  GenomicFeatures::seqlevels(gr) <-  as.character(seqnames(gr))
  gr
})

transcript_coords <- map(gene_ids[1:3], function(x) {
  gr <- transcripts(txdb, filter = list(gene_id = x))
  as.data.frame(gr)
}) %>%
  bind_rows(.id = 'gene')

fls <- list.files('data/bigwig', pattern = 'ercc6_*', full.names = TRUE)

gene_coverage <- map(gene_coords[1:3], function(x) {
  map(fls, function(f) {
    import.bw(con = BigWigFile(f),
              selection = BigWigSelection(x)) %>%
      as.data.frame()
  }) %>%
    bind_rows(.id = 'sample')
}) %>%
  bind_rows(.id = 'gene')

p1 <- list(coverage = gene_coverage %>%
       mutate(x = (start + end)/2) %>%
       filter(score < 2),
     gene_annotation = transcript_coords %>%
       mutate(y = as.numeric(as.factor(tx_name)))) %>%
 bind_rows(.id = 'type') %>%
  ggplot() +
  geom_segment(aes(x = x, y = 0, xend = x + 3, yend = score, group = sample), alpha = .1) +
  geom_segment(aes(x = start, xend = end, y = y, yend = y)) +
  facet_wrap(type~gene, nrow = 2, scales = 'free') +
  theme(legend.position = 'none') +
  theme_void() +
  theme(strip.text = element_blank())

transcript_coords <- map(gene_ids[c(1, 4:6)], function(x) {
  gr <- transcripts(txdb, filter = list(gene_id = x))
  as.data.frame(gr)
}) %>%
  bind_rows(.id = 'gene')

fls <- list.files('data/bigwig', pattern = 'vezf1_*', full.names = TRUE)

gene_coverage <- map(gene_coords[c(1, 4:6)], function(x) {
  map(fls, function(f) {
    import.bw(con = BigWigFile(f),
              selection = BigWigSelection(x)) %>%
      as.data.frame()
  }) %>%
    bind_rows(.id = 'sample')
}) %>%
  bind_rows(.id = 'gene')

p2 <- list(coverage = gene_coverage %>%
        mutate(x = (start + end)/2),
      gene_annotation = transcript_coords %>%
        mutate(y = as.numeric(as.factor(tx_name)))) %>%
    bind_rows(.id = 'type') %>%
    ggplot() +
    geom_segment(aes(x = x, y = 0, xend = x + 3, yend = score, group = sample), alpha = .1) +
    geom_segment(aes(x = start, xend = end, y = y, yend = y)) +
    facet_wrap(type~gene, nrow = 2, scales = 'free') +
    theme(legend.position = 'none') +
    theme_void() +
    theme(strip.text = element_blank())

plot_grid(p1, p2,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10,
          nrow = 2) %>%
  ggsave(plot = .,
         width = 20, height = 24, units = 'cm',
         filename = paste(figures_dir, 'chip_binding.png', sep = '/'))

