#' Get GO annotations
#'
#' Get Gene Ontology (GO) annotations is a tidy data.frame. This function is
#' customized to extract the lower level GO IDs before finding the gene
#' products associated with each ID. Finally, the term corresponding to each
#' ID is added
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
#' @return A data.frame of at least 6 columns
#' \itemize{
#' \item category The category/ies as in go_names
#' \item go The GO ID/s as in go_id or their lower level components
#' \item evidence The type of evidence for each gene product
#' \item ontology The type of ontology (e.g. BP, MF and/or CC)
#' \item (column/s) The query columns as in columns
#' \item term The term corresponding to each GO ID
#' }
#'
#' #' @examples
#' \dontrun{
#' # load required libraries
#' library(org.Mm.eg.db)
#' library(GO.db)
#'
#' # prepare query
#' go_ids <- c('GO:0006914', 'GO:0004679')
#' go_names <- c('autophagy', 'AMPK')
#' go_db <- GOBPCHILDREN
#' go_term <- GOTERM
#' org_db <- org.Mm.eg.db
#' columns <- 'SYMBOL'
#'
#' # call function
#' ann <- annotation_get(go_ids,
#'                       go_names,
#'                       go_db,
#'                       go_term,
#'                       org_db,
#'                       columns)
#' }
#'
#' @importFrom purrr map possibly
#' @importFrom AnnotationDbi mget Term
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

#' Reshape data for WGCNA
#'
#' Reshap an expression matrix to the format required by WGCNA
#'
#' @param mat A numeric expression matrix
#' @param log A logical (default, FALSE) to take the log matrix
#' @param zero_neg A logical (default, FALSE) to zero all negative values
#'
#' @return A list with one data.frame item, data contains the transposed matrix.
#' This is the genes in column and samples in rows
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
#' Get the metadata associated with an ExpressionSet object. This could be
#' the phenotype or the features data. In either cases the functions uses
#' the standard bioconductor accessor with some formating
#'
#' @param eset An ExpressionSet object
#' @param col_name A character vector of columns in the metadata
#' @param new_name A character string of the new names of the columns
#' @param type A character string (default 'phenotype')
#'
#' @return A data.frame always contains a column for each col_name with the
#' name in new_name. When type is 'phenotype', row.names are the samples,
#' and when type is 'features', row.names are the feature IDs/probes
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
#' A wrapper to run WGCNA with certain options in a particular order
#'
#' @param dat A data.frame of transposed expression matrix
#' @param power An integer power of the adjacency matrix
#'
#' @return A list of 9 items:
#' \itemize{
#' \item adj A matrix, the adjacency matrix result of adjacency
#' \item tom A matrix, the TOM similarity matrix result of TOMsimilarity
#' \item diss A matrix, the dissimilary matrix (1 - tom)
#' \item modules A numeric, the assignment of each column to a module after the cutting
#' \item colors A character, the colors given to each integer in modules, result of labels2colors
#' \item mes A data.frame, the eigengenes item of calling moduleEigengenes
#' \item me_dist A matrix, the distance matrix (1 - mes)
#' \item me_tree An hclust, the hierarchial clustering object of me_dist
#' \item merged colors A character, the final lables after mergeCloseModules
#' }
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
#' Awrapper to run comparCluster, gofiler and extract results
#'
#' @param index A list of items each is a set og gene IDs
#' @param level An integer, the Gene Ontology
#' @param ... Other arguments passed to compareCluster
#'
#' @return A data.frame as in the slot compareClusterResult
#'
#' @examples
#' \dontrun{
#' # load required libraries
#' library(org.Mm.eg.db)
#'
#' # make an index of gene sets
#' ind <- list(blue = c("Acbd5", "Atg5", "Atm", "Bmf"),
#'             brown = c("Atg4c", "Il3", "Lamp3", "Lepr", "Mefv"))
#'
#' # run comparison
#' module_compare(ind,
#'                fun = 'enrichGO',
#'                OrgDb = org.Mm.eg.db,
#'                keyType = 'SYMBOL',
#'                ont = 'MF',
#'                pAdjustMethod = 'fdr')
#' }
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
#' A wrapper to get STRING interactions in a nice format
#'
#' @param genes A data.frame of at least one column of gene sybols
#' @param evidence A logical (default FALSE) of whether to include the evidence
#' @param ... Other arguments passed to new
#'
#' @return A data.frame of at least two columns in the form of a edgelist with
#' an additional column, evidence when evidence is TRUE
#'
#' @examples
#' \dontrun{
#' genes = c("Atg4c", "Il3", "Lamp3", "Lepr", "Mefv")
#' interactions_get(genes)
#' }
#'
#' @import STRINGdb
#' @importFrom readr read_delim cols_only
#' @importFrom dplyr inner_join starts_with filter
#' @importFrom tidyr gather
#' @importFrom stats na.omit setNames
#' @importFrom magrittr %>%
#'
#' @export
interactions_get <- function(genes, evidence = FALSE, ...) {

  ptn <- list()
  ptn$new <- STRINGdb$new(...)
  ptn$string_mapped <- ptn$new$map(genes,
                                   'symbol',
                                   removeUnmappedRows = TRUE)
  ptn$string_ids <- ptn$string_mapped$STRING_id
  ptn$interactions <- ptn$new$get_interactions(ptn$string_ids)

  if(evidence == TRUE) {
    df <- ptn$interactions %>%
      inner_join(ptn$string_mapped,
                 by = c(from = "STRING_id")) %>%
      inner_join(ptn$string_mapped,
                 by = c(to = "STRING_id")) %>% na.omit() %>%
      dplyr::select(-from, -to, -starts_with("color")) %>%
      gather(evidence, value, -starts_with("symbol")) %>%
      filter(value != 0) %>%
      setNames(c("from", "to", "evidence", "value"))
  } else {
    df <- ptn$interactions %>%
      dplyr::select(from, to) %>%
      inner_join(ptn$string_mapped,
                 by = c(from = "STRING_id")) %>%
      inner_join(ptn$string_mapped,
                 by = c(to = "STRING_id")) %>%
      na.omit() %>%
      dplyr::select(starts_with("symbol")) %>%
      setNames(c("from", "to"))
  }
  return(df)
}


#' Make graph object of module members
#'
#' Make undirected graphs of a list of modules
#'
#' @param modules A list of module members
#' @param edge_list An edge list of the members connections
#'
#' @return A list of igraph objects
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
#' Calculate several centrality measures on the graphs
#'
#' @param gs A list of graph objects
#' @param measure A character string of the centrality measure
#'
#' @return A data.frame. When measure is 'all'
#' \dontrun{
#' \item module The names in modules
#' \item names The names of the entries in modules
#' \item degree The degree centrality of each node, as in centr_degree
#' \item betweenness The betweennes centrality of each node, as in centr_betw
#' \item closeness The closeness centrality of each node, as in centr_clo
#' \item hub The hub score of each node, as in hub_score
#' }
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
#' @return A data.fame. When measure is 'all'
#' \itemize{
#' \item x1 A gene ID
#' \item x2 A gene ID, as in edge list (x1 --> x2)
#' \item value The value of the measure
#' \item measure either Adjacency, TOM or Pearson
#' }
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
