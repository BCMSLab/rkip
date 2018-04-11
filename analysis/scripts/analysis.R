# load required libraries
library(rkip)
library(tidyverse)
library(WGCNA)
library(reshape2)
library(org.Hs.eg.db)

allowWGCNAThreads(3)

# read data
df <- read_tsv('data/atlas_expression.tsv')
mat <- dplyr::select(df, starts_with('GSM')) %>% as.matrix()
rownames(mat) <- df$DesignElementAccession

# collapse rows
mat <- collapseRows(mat, rowID = rownames(mat), rowGroup = df$`Gene Name`)[[1]]

# get phenotype data
design <- read_tsv('data/atlas_design.tsv') %>%
  filter(Assay %in% colnames(mat)) %>%
  dplyr::select(Assay, `Sample Characteristic[disease]`) %>%
  setNames(c('sample', 'disease'))

all(colnames(mat) %in% design$sample)

# get annotatios
ann <- read_tsv('data/annotations.tsv')
ind <- rownames(mat) %in% unique(ann$symbol)

# subset matrix
dat <- t(mat[ind,])

# pick threshold
sft <- pickSoftThreshold(dat)

# run cna
net <- cna_run(dat, power = 5)

# get module members
genes <- data_frame(symbol = colnames(dat),
                    color = net$merged_colors)

# preservation
## read other studies
fls <- list.files('data', pattern = 'prad_*', full.names = TRUE)
names(fls) <- str_split(fls, '/|\\.', simplify = TRUE)[, 2]
md <- map(fls, function(x) {
  df <- read_tsv(x)
  ind <- goodSamplesGenes(df, minFraction = 0)
  df <- df[ind$goodSamples, ind$goodGenes]
  list(data = df)
})
md$ref <- list(data = dat)

mc <- net$merged_colors
names(mc) <- names(dat)
mc <- list(ref = mc)

module_preserv <- modulePreservation(multiData = md,
                                     multiColor = mc,
                                     referenceNetworks = 9,
                                     savePermutedStatistics = FALSE,
                                     nPermutations = 5)

# clean and save
rm(df, ind, mc, fls)

save.image('data/rkip_wgcna.rda')
