#' Sample to sample distances
#'
#' Plot sample distances between barcode sets / samples
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param method Distance metric to use (string). One of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Refer to stats::dist for details. Default = `euclidean`.
#' @param upper Plot only the upper half of the matrix (boolean). Default = `TRUE`.
#' @param clustered Cluster rows and columns of matrix (boolean). Default = `TRUE`.
#' @param title Title of plot (string). Default = `Sample Distances - method`.
#'
#' @return Returns a heatmap of sample distances using desired clustering
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' plotBarcodeDistance(test.dge)

plotBarcodeDistance <- function(dgeObject,
                                method = "euclidean",
                                upper = TRUE,
                                clustered = TRUE,
                                title = "Sample Distances") {
  if (methods::is(dgeObject)[1] != "DGEList") {
    stop("Please supply a valid DGEList object as input")
  }
  counts <- dgeObject$counts

  # generate correlation matrix
  sample.distances <- stats::dist(t(counts), method = method)
  sample.matrix <- as.matrix(sample.distances)

  # cluster matrix if specified
  if (isTRUE(clustered)) {
    sample.matrix <- cluster_cormat(sample.matrix)
  }

  # take only upper triangle if specified
  if (isTRUE(upper)) {
    sample.matrix <- get_upper_tri(sample.matrix)
  }

  # Melt the correlation matrix
  melted_dist <- reshape2::melt(sample.matrix, na.rm = TRUE)

  # Create a ggheatmap
  title <- paste(title, "-", method)
  ggheatmap <-
    ggplot2::ggplot(melted_dist, ggplot2::aes(`Var2`, `Var1`, fill = `value`)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(option = "magma") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 1,
        size = 6,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 6)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "", y = "", title = title)

  return(ggheatmap)
}
