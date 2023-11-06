#' Sample to sample correlation
#'
#' Plot sample correlation between barcode sets / samples
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param method Correlation metric to use (string). One of "pearson", "kendall", or "spearman". Refer to stats:cor for details. Default = `pearson`.
#' @param upper Plot only the upper half of the matrix (boolean). Default = `TRUE`.
#' @param clustered Cluster rows and columns of matrix (boolean). Default = `TRUE`.
#' @param title Title of plot (string). Default = `Sample Correlation - method`.
#'
#' @return Returns a heatmap of sample distances using desired clustering
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' plotBarcodeCorrelation(test.dge)

plotBarcodeCorrelation <- function(dgeObject,
                                   method = "pearson",
                                   upper = TRUE,
                                   clustered = TRUE,
                                   title = "Sample Correlation") {
  inputChecks(dgeObject)
  counts <- dgeObject$counts

  # generate correlation matrix
  cormat <- round(stats::cor(counts, method = method), 2)

  # cluster matrix if specified
  if (isTRUE(clustered)) {
    cormat <- cluster_cormat(cormat)
  }

  # take only upper triangle if specified
  if (isTRUE(upper)) {
    cormat <- get_upper_tri(cormat)
  }

  # Melt the correlation matrix
  melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)

  # Create a ggheatmap
  title <- paste(title, "-", method)
  ggheatmap <-
    ggplot2::ggplot(melted_cormat, ggplot2::aes(`Var2`, `Var1`, fill = `value`)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c() +
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

  # Print the heatmap
  print(ggheatmap)
}
