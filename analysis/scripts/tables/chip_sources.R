# load libraries
library(tidyverse)
library(xtable)

# global variables
tables_dir <- 'manuscript/tables'

'http://cistrome.org/ap/u/mahmoudahmed/h/rkipcommonregulators'

# generate table
fm <- data_frame(file_name = c('49027.bw', '49036.bw', '49033.bw', '49030.bw', '44357.bw', '44356.bw', '8533.bw', '8532.bw', '38706.bw'),
           pmid = c(25249633, 25249633, 25249633, 25249633, 26578602, 26578602, 22308494, 22308494, 23953112),
           geo = c('GSM1405140', 'GSM1405140', 'GSM1405137', 'GSM1405134', 'GSM1215093', 'GSM1215092', 'GSM866197', 'GSM866196', 'GSM1208807'))

sm <- data_frame(pmid = c(25249633, 26578602, 22308494, 23953112),
                 cell_type = c('Fibroblast', 'Fibroblast', 'HeLa','LoVo'),
                 factor = c('ERCC6','ERCC6','VEZF1','VEZF1'))

pm <- data_frame(pmid = c(25249633, 26578602, 22308494, 23953112),
                 citation = c('Wang2014DysregulationDisease.', 'Lake2016TheStress.', 'Gowher2012Vezf1II.', 'Yan2013TranscriptionSites'))

sm %>%
  full_join(fm) %>%
  full_join(pm) %>%
  group_by(factor, cell_type, citation) %>%
  summarise(samples = paste(unique(geo), collapse = ', ')) %>%
  ungroup() %>%
  mutate(citation = paste0('\\cite{', citation, '}')) %>%
  dplyr::select(factor, cell_type, samples, citation) %>%
  setNames(c('Factor', 'Cell Type', 'GEO Accession', 'Reference')) %>%
  xtable(caption = '\\hl{Transcription factors ChIP-Seq data sources.}',
         label = 'tab:chip_sources',
         align = 'clllc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        file = paste(tables_dir, 'chip_sources.tex', sep = '/'))

