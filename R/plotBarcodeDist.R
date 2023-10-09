#' Sample to sample distances
#'
#' Plot sample distances between barcode sets / samples
#'
#' @param counts matrix/dataframe containing raw or normalised counts
#' @param method distance metric to use. Refer to stats::dist for available options
#' @param upper Logical. plot only the upper half of the matrix
#' @param clustered Logical. cluster rows and columns of matrix
#' @param name title of plot
#'
#' @return Returns a heatmap of sample distances using desired clustering
#'
#' @export
#'
#' @examples
#' plotBarcodeDistance(test.dge$counts)

plotBarcodeDistance <- function(counts,
                          method = "euclidean",
                          upper = T,
                          clustered = T,
                          name = "Sample Distances"){

  # generate correlation matrix
  sample.distances <- stats::dist(t(counts), method = method)
  sample.matrix <- as.matrix(sample.distances)

  # cluster matrix if specified
  if(isTRUE(clustered)){
    sample.matrix <- cluster_cormat(sample.matrix)
  }

  # take only upper triangle if specified
  if(isTRUE(upper)){
    sample.matrix <- get_upper_tri(sample.matrix)
  }

  # Melt the correlation matrix
  melted_dist <- reshape2::melt(sample.matrix, na.rm = TRUE)

  # Create a ggheatmap
  name <- paste(name, "-", method)
  ggheatmap <- ggplot2::ggplot(melted_dist, ggplot2::aes(`Var2`, `Var1`, fill = `value`)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(option = "magma") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 6, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "", y = "", title = name)

  # Print the heatmap
  return(ggheatmap)
}

#' Sample to sample correlation
#'
#' Plot sample correlation between barcode sets / samples
#'
#' @param counts matrix/dataframe containing raw or normalised counts
#' @param method correlation metric to use. Refer to stats::cor for available options
#' @param upper Logical. plot only the upper half of the matrix
#' @param clustered Logical. cluster rows and columns of matrix
#' @param name title of plot
#' @param samples Dataframe containing sample metadata sheet with samples as row names.
#' @param group metadata field to annotate samples
#' @param col_annot list of vectors assigning colours to each level of metadata
#'
#' @return Returns a heatmap of sample distances using desired clustering
#'
#' @export
#'
#' @examples
#' test.mat <- matrix(rnorm(10*10,mean=0,sd=2), 10, 10)
#' plotBarcodeCorrelation(counts = test.mat, method = "pearson")

plotBarcodeCorrelation <- function(counts,
                            method = "pearson",
                            upper = T,
                            clustered = T,
                            group = NULL,
                            col_annot = NULL,
                            samples = NULL,
                            name = "Sample Correlation"){

  # check input parameters
  if (!is.null(group)) {
    if (any(is.null(samples), is.null(col_annot))) {
      stop("When specifying a group, please also provide sample metadata and color annotations.")
    }
  }

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

  cormat <- cormat[rev(rownames(cormat)), ]

  if (!is.null(group)) {
    # keep group order for legend in plot
    # reorder samples in case they were clustered
    group_cols <- group[group %in% colnames(samples)]
    df_annotation <-
      samples[colnames(cormat), group_cols, drop = FALSE]

    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = df_annotation,
      which = 'col',
      col = col_annot,
      annotation_width = unit(c(1, 4), 'cm'),
      gap = unit(1, 'mm')
    )
  } else {
    ha <- NULL
  }

  ComplexHeatmap::Heatmap(
    cormat,
    cluster_columns = F,
    cluster_rows = F,
    col = viridis::viridis(100),
    bottom_annotation = ha,
    width = ncol(cormat) * unit(5, "mm"),
    height = nrow(cormat) * unit(5, "mm"),
    heatmap_legend_param = list(title = paste(method, "\ncorrelation")),
    column_title = name,
    na_col = "white",
    row_names_side = "left"
  )
}


#' Cluster correlation matrix
#'
#' cluster a correlation matrix using hierarchical clustering
#'
#' @param cormat matrix of correlation values
#'
#' @return Returns a matrix of correlation values with columns and rows hierarchically clustered
#'
#' @export
#'
cluster_cormat <- function(cormat){
  dd <- stats::as.dist((1-cormat)/2)
  hc <- stats::hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}

#' Get upper triangle
#'
#' change lower triangle values in a mirrored matrix to NA
#'
#' @param cormat matrix of correlation values
#'
#' @return Returns a matrix of correlation values with lower triangle values changed to NA
#' @export
#'
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}
