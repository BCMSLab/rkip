# load libraries
library(tidyverse)
library(rkip)
library(org.Hs.eg.db)
library(cgdsr)

# global variable

# download data
## atlas expression data
download.file(
  url = 'https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-3325/resources/DifferentialSecondaryDataFiles.Microarray/normalized-expressions',
  destfile = 'data/atlas_expression.tsv'
)

download.file(
  url = 'https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-3325/resources/ExperimentDesignFile.Microarray/experiment-design',
  destfile = 'data/atlas_design.tsv'
)


## annotation data
ann <- annotation_get(
  go_id = c('GO:0001837', 'GO:0006914', 'GO:0008429'),
  go_names = c('emt', 'autophagy', 'pebp'),
  org_db = org.Hs.eg.db,
  columns = 'SYMBOL',
  remove_predicted = FALSE
)

write_tsv(ann, 'data/annotations.tsv')

# get cancer studies
cgd <- CGDS('http://www.cbioportal.org/')

#studies_ids <- getCancerStudies(cgd) %>%
#  filter(grepl('prad', cancer_study_id)) %>%
#  pull(cancer_study_id)

#map(studies_ids,
#            function(x) {
#              getGeneticProfiles(cgd,
#                                 cancerStudy = x)
#            })
#map(studies_ids,
#         function(x) {
#           getCaseLists(cgd,
#                        cancerStudy = x)
#         })

study_profile <- list(
  prad_su2c_2015 = 'prad_su2c_2015_rna_seq_mrna',
  prad_broad_2013 = 'prad_broad_2013_mrna',
  prad_broad = 'prad_broad_mrna',
  prad_fhcrc = 'prad_fhcrc_rna_seq_mrna',
  prad_mskcc = 'prad_mskcc_mrna_zbynorm',
  prad_tcga_pub = 'prad_tcga_pub_rna_seq_v2_mrna',
  prad_tcga = 'prad_tcga_rna_seq_v2_mrna',
  prad_mskcc_cheny1_organoids_2014 = 'prad_mskcc_cheny1_organoids_2014_rna_seq_rna'
)

gene_ids <- unique(ann$symbol)

imap(study_profile, function(x, .y) {
  getProfileData(cgd,
                 genes = gene_ids,
                 geneticProfiles = x,
                 caseList = paste(.y, 'all', sep = '_')) %>%
    write_tsv(path = paste('data/', .y, '.tsv', sep = ''))
  return(NULL)
})

# get string interacitons
df <- interactions_get(data.frame(symbol = unique(ann$symbol)),
                       input_directory = 'data/',
                       species = 9606,
                       evidence = TRUE)
write_tsv(df, 'data/string_interactions.tsv')
