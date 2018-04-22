# load libraries
library(tidyverse)
library(rkip)
library(org.Hs.eg.db)
library(cgdsr)
library(cRegulome)

# download data
## atlas expression data
if(!file.exists('data/atlas_expression.tsv')) {
  download.file(
    url = 'https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-3325/resources/DifferentialSecondaryDataFiles.Microarray/normalized-expressions',
    destfile = 'data/atlas_expression.tsv'
  )
}

if(!file.exists('data/atlas_design.tsv')) {
  download.file(
    url = 'https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-3325/resources/ExperimentDesignFile.Microarray/experiment-design',
    destfile = 'data/atlas_design.tsv'
  )
}

## annotation data
if(!file.exists('data/annotations.tsv')) {
  ann <- annotation_get(
    go_id = c('GO:0001837', 'GO:0006914', 'GO:0008429'),
    go_names = c('emt', 'autophagy', 'pebp'),
    org_db = org.Hs.eg.db,
    columns = 'SYMBOL',
    remove_predicted = FALSE
  )

  write_tsv(ann, 'data/annotations.tsv')
}

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

ann <- read_tsv('data/annotations.tsv')
gene_ids <- unique(ann$symbol)

imap(study_profile, function(x, .y) {
  fl <- paste('data/', .y, '.tsv', sep = '')
  if(!file.exists(fl)) {
    getProfileData(cgd,
                   genes = gene_ids,
                   geneticProfiles = x,
                   caseList = paste(.y, 'all', sep = '_')) %>%
      write_tsv(path = fl)
  }
  return(NULL)
})

# get string interacitons
if(!file.exists('data/string_interactions.tsv')) {
  df <- interactions_get(data.frame(symbol = unique(ann$symbol)),
                         input_directory = 'data/',
                         species = 9606,
                         evidence = TRUE)
  write_tsv(df, 'data/string_interactions.tsv')
}

# get tf data
#tf <- c('ERCC6', 'VEZF1')
#map(tf, function(x) {
#  fl <- paste('data/', x, '_cor.tsv', sep = '')
#  if(!file.exists(fl)) {
#    url <- paste('http://cistrome.org/CistromeCancer/CancerTarget/examples/',
#                 x, '.cor.csv',
#                 sep = '')
#    read_csv(url) %>%
#      write_tsv(path = fl)
#  }
#  return(NULL)
#})

#map(tf, function(x) {
#  fl <- paste('data/', x, '_rp.tsv', sep = '')
#  if(!file.exists(fl)) {
#    url <- paste('http://cistrome.org/CistromeCancer/CancerTarget/examples/',
#                 x, '.rp.csv',
#                 sep = '')
#    read_csv(url) %>%
#      write_tsv(path = fl)
#  }
#  return(NULL)
#})

# predicted miRNA
#url <- 'http://www.targetscan.org/vert_72/temp/TargetScan_7.2_ENST00000261313.2_predicted_targeting_details.txt'
#if(!file.exists('data/predicted_targeting.tsv')) {
# download.file(url,
#               destfile = 'data/predicted_targeting.tsv')
#}

if(!file.exists('data/cRegulome.db')) {
  get_db(test = FALSE)
}
