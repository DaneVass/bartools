#' @title
#' Group barcodes to clones
#'
#' @description
#' Group barcodes to clones by connected components given a similarity matrix and threshold.
#'
#' @param bcSimilarity Similarity matrix with barcodes as rows and columns (output of calculateBarcodeSimilarities).
#' @param threshold Similarity metric threshold above which barcodes will be grouped together (decimal). Default = `0.7`.
#'
#' @return A named list with barcodes as names and clone identity as values.
#' @export
#'

# group barcodes into clones by connected components and a similarity threshold
barcodesToClones <- function(bcSimilarity, threshold = 0.7) {
  print(paste0(
    sum(bcSimilarity >= threshold, na.rm = T),
    " barcodes will be collapsed into clones"
  ))
  bc_similarity_binary <- bcSimilarity >= threshold
  bc_similarity_graph <-
    igraph::graph_from_adjacency_matrix(bc_similarity_binary)
  clones <- igraph::components(bc_similarity_graph)$membership
  return(clones)
}
