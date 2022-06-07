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
#' test.mat <- matrix(rnorm(10*10,mean=0,sd=2), 10, 10)
#' plotBarcodeDistance(test.mat)

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
  ggheatmap <- ggplot2::ggplot(melted_dist, ggplot2::aes(melted_dist$Var2, melted_dist$Var1, fill = melted_dist$value)) +
    ggplot2::geom_tile(color = "white") +
    viridis::scale_fill_viridis(option = "inferno") +
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
                            name = "Sample Correlation"){

  # generate correlation matrix
  cormat <- round(stats::cor(counts, method = method),2)

  # cluster matrix if specified
  if(isTRUE(clustered)){
    cormat <- cluster_cormat(cormat)
  }

  # take only upper triangle if specified
  if(isTRUE(upper)){
    cormat <- get_upper_tri(cormat)
  }

  # Melt the correlation matrix
  melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)

  # Create a ggheatmap
  name <- paste(name, "-", method)
  ggheatmap <- ggplot2::ggplot(melted_cormat, ggplot2::aes(melted_cormat$Var2, melted_cormat$Var1, fill = melted_cormat$value)) +
    ggplot2::geom_tile(color = "white") +
    viridis::scale_fill_viridis() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 6, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "", y = "", title = name)

  # Print the heatmap
  print(ggheatmap)
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
