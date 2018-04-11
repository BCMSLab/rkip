#' Get GO annotations
#'
#' @param go_id A character vector of GO IDs
#' @param go_names A character vector of GO category names
#' @param go_db A GO.db object
#' @param go_term A GO.db object for GO terms
#' @param org_db An organism object
#' @param columns A charachter vector
#' @param remove_predicted A logical (default, TRUE) for removing predicted genes
#' @param ... Other arguments passed to select
#'
#' @return A data.frame
#'
#' @importFrom purrr map possibly
#' @importFrom AnnotationDbi select mget Term
#' @importFrom dplyr bind_rows filter mutate
#' @importFrom magrittr %>%
#'
#' @export
annotation_get <- function(go_id, go_names, go_db, go_term, org_db, columns, remove_predicted = TRUE, ...) {
  # make a list of ids
  ids <- (go_id)

  # get ids of certain level when applicable
  if(!missing(go_db)) {
    ids <- ids %>%
      map(possibly(function(x) mget(x, go_db) %>% unlist, NA))

    ind <- is.na(ids)
    ids[ind] <- go_id[ind]
  }

  names(ids) <- go_names

  # get annotation columns
  df <- map(ids, function(x) {
    AnnotationDbi::select(org_db,
                          x,
                          columns,
                          'GO')
  }) %>%
    bind_rows(.id = 'category') %>%
    na.omit()

  if(!missing(go_term)) {
    df <- df %>%
      mutate(term = sapply(mget(GO, envir = go_term),
                           FUN = Term))
  }

  names(df) <- tolower(names(df))
  if(remove_predicted == TRUE) {
    df <- df %>%
      filter(!grepl('Gm*', symbol)) %>%
      na.omit()
  }
  return(df)
}

#' Subset an ExpressionSet
#'
#' @param eset An ExpressionSet object
#' @param byrow A logical vector masking rows of eset
#' @param bycol A logical vector masking columns of eset
#' @param collapse_rows A logical (default FALSE)
#' @param ... Other arguments passed to collapse_rows
#'
#' @return A tidy expression matrix
#'
#' @examples
#' library(Biobase)
#' data("sample.ExpressionSet")
#'
#' pd <- pData(sample.ExpressionSet)
#'
#' set.seed(123)
#' byrow <- sample(c(TRUE, FALSE), nrow(sample.ExpressionSet), replace = TRUE)
#' bycol <- pd$sex == 'Female'
#'
#' expression_subset(sample.ExpressionSet, byrow, bycol)
#'
#' @importFrom Biobase exprs
#' @importFrom WGCNA collapseRows
#'
#' @export
expression_subset <- function(eset, byrow, bycol, collapse_rows=FALSE, ...) {

  if(missing(byrow)) byrow <- TRUE
  if(missing(bycol)) bycol <- TRUE

  if(nrow(eset) != length(byrow) & byrow != TRUE) stop('byrow needs to be same length as nrow(eset).')
  if(ncol(eset) != length(bycol) & bycol != TRUE) stop('bycol needs to be same length as ncol(eset).')

  e <- eset[byrow, bycol]
  mat <- exprs(e)

  if(collapse_rows == TRUE) {
    #if(missing(rowGroup)) stop('rowGroup is required for collapse Rows.')
    #if(missing(rowID)) stop('rowID is required for collapse Rows.')

    #if(length(rowGroup) != length(rowID)) stop('rowGroup and rowID need to be of same length.')

    collapsed_mat <- collapseRows(mat, ...)
    mat <- collapsed_mat[[1]]
  }
  if(sum(duplicated(rownames(mat)))) warning('matrix rows are not unique')

  return(mat)
}

#' Reshape data for WGCNA
#'
#' @param mat A numeric expression matrix
#' @param log A logical (default, FALSE) to take the log matrix
#' @param zero_neg A logical (default, FALSE) to zero all negative values
#'
#' @return A list with one data.frame slot, data contains the transposed matrix
#'
#' @examples
#' library(Biobase)
#' data(sample.ExpressionSet)
#'
#' mat <- exprs(sample.ExpressionSet)
#'
#' expression_reshape(mat)
#'
#' @export
expression_reshape <- function(mat, log = FALSE, zero_neg = FALSE) {
  if(zero_neg == TRUE) {
    mat[mat < 0] <- 0
  }

  if(log == TRUE) {
    mat <- log(mat + 1)
  }

  # return
  list(data = t(mat))
}

#' Get metadata from ExpressionSet
#'
#' @param eset An ExpressionSet object
#' @param col_name A character vector of columns in the metadata
#' @param new_name A character string of the new names of the columns
#' @param type A character string (default 'phenotype')
#'
#' @return A data.frame
#'
#' @examples
#' library(Biobase)
#' data("sample.ExpressionSet")
#'
#' metadata_get(sample.ExpressionSet,
#'              c('sex', 'type'),
#'              c('sex', 'type'))
#'
#' @importFrom Biobase pData fData sampleNames featureNames
#'
#' @export
metadata_get <- function(eset, col_name, new_name, type = 'phenotype') {
  switch (type,
          'phenotype' = {
            df <- as.data.frame(pData(eset)[, col_name])
            if(!missing(eset)) {
              names(df) <- new_name
            }
            df$gsm <- sampleNames(eset)
            return(df)
          },
          'features' = {
            df <- as.data.frame(fData(eset)[, col_name])
            if(!missing(new_name)) {
              names(df) <- new_name
            }
            df$probe_id <- featureNames(eset)
            return(df)
          }
  )
}

#' Run WGCNA
#'
#' @param dat A data.frame of transposed expression matrix
#' @param power An integer power of the adjacency matrix
#'
#' @return A list
#'
#' @examples
#' library(Biobase)
#' data(sample.ExpressionSet)
#' dat <- list(data = t(exprs(sample.ExpressionSet)))
#'
#' cna_run(dat$data, power = 5)
#'
#' @import WGCNA
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats as.dist hclust
#'
#' @export
cna_run <- function(dat, power) {

  # make an empty list
  net <- list()

  # adjacency matrix
  net$adj <- adjacency(dat, power = power)

  # tom similarity
  net$tom <- TOMsimilarity(net$adj)

  # dissimilarity matrix
  net$diss <- 1 - net$tom

  # cut modules
  net$modules <- cutreeDynamic(hclust(as.dist(net$diss), method = 'average'),
                               distM = net$diss,
                               deepSplit = FALSE,
                               pamRespectsDendro = FALSE,
                               minClusterSize = 10)

  # labels to colors
  net$colors <- labels2colors(net$modules)

  # calculate module eigengenes
  net$mes <- moduleEigengenes(dat, colors = net$colors)$eigengenes

  # merge module
  net$me_dist <- 1-cor(net$mes)
  net$me_tree <- hclust(as.dist(net$me_dist), method = 'average')

  net$merged_colors <- mergeCloseModules(dat,
                                         net$colors,
                                         cutHeight = .6)$colors
  net$mes <- orderMEs(moduleEigengenes(dat, net$merged_colors)$eigengenes)

  return(net)
}

#' Compare enrichment by modules
#'
#' @param index A list
#' @param level An integer
#' @param ... Other arguments passed to compareCluster
#'
#' @return A data.frame
#'
#' @importFrom clusterProfiler compareCluster gofilter
#'
#' @export
module_compare <- function(index, level = 4, ...) {
  comp <- compareCluster(index, ...)
  comp <- gofilter(comp, level)
  comp@compareClusterResult
}
#' Get STRING interactions
#'
#' @param genes A data.frame of at least one column of gene sybols
#' @param evidence A logical (default FALSE) of whether to include the evidence
#' @param ... Other arguments passed to new
#'
#' @return A data.frame
#'
#' @import STRINGdb
#' @importFrom readr read_delim cols_only
#' @importFrom dplyr select inner_join starts_with filter
#' @importFrom tidyr gather
#' @importFrom stats na.omit setNames
#' @importFrom magrittr %>%
#'
#' @export
interactions_get <- function(genes, input_directory, evidence = FALSE, ...) {

  fl <- list.files(input_directory, pattern = 'protein_links.tsv.gz', full.names = TRUE)

  ptn <- list()
  ptn$new <- STRINGdb$new(input_directory = input_directory, ...)
  ptn$string_mapped <- ptn$new$map(genes, 'symbol', removeUnmappedRows = TRUE)
  ptn$string_ids <- ptn$string_mapped$STRING_id

  if(!file.exists(fl)) {
    ptn$new$get_interactions(ptn$string_ids)
  }

  if(evidence == TRUE) {
    ptn$interactions <- read_delim(fl, delim = ' ')

    df <- ptn$interactions %>%
      inner_join(ptn$string_mapped, by = c('protein1' = 'STRING_id')) %>%
      inner_join(ptn$string_mapped, by = c('protein2' = 'STRING_id')) %>%
      na.omit() %>%
      select(-starts_with('protein'), -starts_with('color')) %>%
      gather(evidence, value, -starts_with('symbol')) %>%
      filter(value != 0) %>%
      setNames(c('from', 'to', 'evidence', 'value'))
  } else {
    ptn$interactions <- read_delim(fl, delim = ' ',
                                   col_types = cols_only(protein1 = 'c', protein2 = 'c'))
    df <- ptn$interactions %>%
      inner_join(ptn$string_mapped, by = c('protein1' = 'STRING_id')) %>%
      inner_join(ptn$string_mapped, by = c('protein2' = 'STRING_id')) %>%
      na.omit() %>%
      select(starts_with('symbol')) %>%
      setNames(c('from', 'to'))
  }
  return(df)
}


#' Make graph object of module members
#'
#' @param modules A list of module members
#' @param edge_list An edge list of the members connections
#'
#' @return A list of graph objects
#'
#' @importFrom purrr map
#' @importFrom igraph graph_from_data_frame subgraph V
#'
#' @export
module_network <- function(modules, edge_list) {
  map(modules, function(x) {
    edge_list %>%
      graph_from_data_frame(directed = FALSE) %>%
      subgraph(v = V(.)$name %in% x)
  })
}

#' Calculate member importance in a graph
#'
#' @param gs A list of graph objects
#' @param measure A character string of the centrality measure
#'
#' @return A data.frame
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @import igraph
#'
#' @export
member_importance <- function(gs, measure = 'all') {
  map(gs, function(g) {
    list(names = get.vertex.attribute(g, 'name'),
         degree = centr_degree(g, normalized = TRUE)$res,
         betweenness = centr_betw(g, directed = FALSE, normalized = TRUE)$res,
         closeness = centr_clo(g, normalized = TRUE)$res,
         hub = hub_score(g)$vector) %>%
      bind_rows()
  }) %>%
    bind_rows(.id = 'module')
}

#' Calculate pair-wise similarity in an expresison matrix
#'
#' @param dat A data.frame of transposed expression matrix
#' @param measure A character string of the similarity measure
#'
#' @return A data.fame
#'
#' @importFrom WGCNA adjacency TOMsimilarity cor
#' @importFrom reshape2 melt
#' @importFrom dplyr filter
#' @importFrom stats setNames
#'
#' @export
member_similarity <- function(dat, measure = 'all') {
  sim <- list()
  sim$Adjacency <- adjacency(dat, power = 5)
  sim$TOM <- TOMsimilarity(sim$Adjacency)
  attr(sim$TOM, 'dimnames') <- attr(sim$Adjacency, 'dimnames')
  sim$Pearson <- cor(dat)

  melt(sim) %>%
    setNames(c('x1', 'x2', 'value', 'measure'))
}
