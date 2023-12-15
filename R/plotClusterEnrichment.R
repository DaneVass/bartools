#' plotClusterEnrichment
#'
#' Hypergeometric test for enrichment of specific classes of cells in
#' single cell RNA-seq clusters
#'
#' @title
#' Hypergeometric test for cluster enrichment, e.g. to test whether a barcode is enriched in any cluster.
#'
#' @description
#' Takes a single cell object, a grouping variable, a factor within the group to test,
#' and an ident class (i.e. clusters). Per level of ident performs hypergeometric testing
#' for enrichment of group factor
#'
#' @param sc.obj single cell object in Seurat or SingleCellExperiment format containing clusters and group metadata
#' @param group a column of metadata (string). Default = `barcode`.
#' @param factor a level of group to test for enrichment per cluster, e.g. a specific barcode (string).
#' @param clusters a column of metadata defining the cluster identities of each cell (string). Default = `seurat_clusters`
#' @param threshold P-value threshold for hypergeometric test (decimal). Default = `0.01`.
#' @param order Order clusters in plot by -log10 P-value (boolean). Default = `TRUE`.
#' @param plot Plot hypergeometric test results or return plot data (boolean). Default = `TRUE`.
#'
#' @return Returns a histogram of hypergeometric test results or underlying plot data
#' @export
#'

plotClusterEnrichment <- function(sc.obj = NULL,
                                  group = "barcode",
                                  factor = NULL,
                                  clusters = "seurat_clusters",
                                  threshold = 0.01,
                                  order = TRUE,
                                  plot = TRUE) {

  # check inputs
  if (is.null(sc.obj)) {
    stop("Please supply a Seurat or SingleCellExperiment object")
  }

  # get metadata and ident class
  if (class(sc.obj)[1] == "Seurat") {
    meta <- sc.obj@meta.data
    type <- "Seurat"
  } else {
    if (class(sc.obj)[1] == "SingleCellExperiment") {
      meta <- as.data.frame(sc.obj@colData)
      type <- "SingleCellExperiment"
    } else {
      stop("Single cell object must be supplied in Seurat or SingleCellExperiment format")
    }
  }

  # input checks

  # group must be column in metadata
  if (!group %in% colnames(meta)) {
    stop("Group variable is not column in metadata")
  }
  # factor must be in group column
  if (!factor %in% meta[group]) {
    stop("Factor is not present in group column")
  }
  # clusters must be column in metadata
  if (!clusters %in% colnames(meta)) {
    stop("Clusters variable is not column in metadata")
  }
  if (length(factor) > 1) {
    stop("Can only test one factor at a time")
  }

  # get total cells and test cells
  total.cells <- ncol(sc.obj)

  factor.rows <- which(meta[,`group`, drop = TRUE] == as.character(factor))
  test.cells <- rownames(meta[factor.rows,])
  test.cells.total <- length(test.cells)
  clusters <- levels(meta[,`clusters`])

  # loop phyper tests for factor enrichment over each level of clusters
  cluster.hyper.tests <- lapply(clusters, function(x) {
    cluster.all.rows <- which(meta[,`clusters`, drop = TRUE] == as.character(x))
    cluster.total.cells <- nrow(meta[cluster.all.rows,])
    cluster.test.rows <- which(meta[cluster.all.rows,`group`, drop = TRUE] == as.character(factor))
    cluster.test.cells <- nrow(meta[cluster.test.rows,])
    message("---")
    message(paste("cluster", x, sep = "_"))
    message(paste("all cells:", total.cells, sep = " "))
    message(paste("all", as.character(factor), "cells:", test.cells.total, sep = " "))
    message(paste("universe:", total.cells - test.cells.total, sep = " "))
    message(paste("cluster total cells:", cluster.total.cells, sep = " "))
    message(paste("cluster", as.character(factor), "cells:", cluster.test.cells, sep = " "))

    # run hypergeometric test for the cluster
    p <- stats::phyper(cluster.test.cells,
                       test.cells.total,
                       total.cells - test.cells.total,
                       cluster.total.cells,
                       lower.tail = F, log.p = FALSE)

    print(paste("Hypergeometric test p-value:", p, sep = " "))
    return(p)
  }
  )

  # generate output data
  names(cluster.hyper.tests) <- paste("cluster", clusters, sep = "_")
  plot.dat <- data.frame(unlist(cluster.hyper.tests),
                         row.names = names(cluster.hyper.tests),
                         neglog10pval = -log10(unlist(cluster.hyper.tests)),
                         enriched = ifelse(unlist(cluster.hyper.tests) < threshold, "YES", "NO"))

  # reorder rows based on pvalue
  if (order) {
    plot.dat <- dplyr::mutate(plot.dat,
                              cluster = factor(rownames(plot.dat),
                                               levels = rownames(plot.dat[order(plot.dat$neglog10pval, decreasing = F),])))
  } else {
    plot.dat <- dplyr::mutate(plot.dat, cluster = factor(rownames(plot.dat)))

  }

  # plot p value histogram or return raw data
  if (plot) {
    p <- ggplot2::ggplot(plot.dat,
                         ggplot2::aes(y = cluster, x = neglog10pval, fill = enriched)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_vline(xintercept = -log10(threshold)) +
      ggplot2::scale_x_continuous(trans = "log1p") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_fill_manual(values = c("blue2", "red2")) +
      ggplot2::ggtitle(paste("Hypergeometric test for", as.character(factor), "enrichment per cluster"))

    return(p)

  } else {
    return(plot.dat)
  }
}
