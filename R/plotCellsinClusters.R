## plotCellsInClusters
#'
#' Plot number or percentage of a group per level of ident in a single
#' cell object
#'
#' @title
#' Plot distribution of group across clusters
#'
#' @description
#' Takes a single cell object, a grouping variable, a factor within the group to test,
#' and an ident class (i.e. clusters). Plots percentage or raw number of cells within factor
#' present per level of ident (i.e. per cluster)
#'
#' @param sc.obj Single cell Seurat or SingleCellExperiment object containing clusters and group metadata.
#' @param group a column of metadata (string). Default = `barcode`.
#' @param factor a level of group to test for enrichment per cluster, e.g. a specific barcode (string).
#' @param clusters a column of metadata defining the cluster identities of each cell (string). Default = `seurat_clusters`
#' @param plot.pct Plot percentages (`TRUE`) or raw numbers (`FALSE`) (boolean). Default = `TRUE`
#' @param plot Create plot (`TRUE`) or return data (`FALSE`) (boolean). Default = `TRUE`.
#'
#' @importFrom rlang .data
#' @return Returns a histogram or underlying plot data
#' @export


plotCellsInClusters <- function(sc.obj,
                                group = "barcode",
                                factor = NULL,
                                clusters = "seurat_clusters",
                                plot.pct = T,
                                plot = T) {

  idents <- clusters
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
      meta <- sc.obj@colData
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
  if (!factor %in% meta[[group]]) {
    stop("Factor is not present in group column")
  }
  # clusters must be column in metadata
  if (!idents %in% colnames(meta)) {
    stop("Clusters variable is not column in metadata")
  }
  if (length(factor) > 1) {
    stop("Can only test one factor at a time")
  }

  # get number of seurat clusters. NB plotting is 0 indexed so total - 1
  clusters <-
    data.frame(clusters = as.factor(seq(0, max(
      as.numeric(meta[, idents]) - 1
    ), 1)))
  colnames(clusters) <- as.character(idents)
  # get and count cells of interest
  meta.bc <- as.data.frame(meta[!is.na(meta[, `group`]), ])
  # get cells in group that are factor
  meta.select <- meta.bc[which(meta.bc[, `group`] == factor), ]

  # group by and tally clusters (i,e. idents)
  idents <- sym(idents)
  meta.bc.group <-
    dplyr::tally(dplyr::group_by(meta.select,!!idents))
  print(meta.bc.group)

  # join counts per cluster (ident) of group factor into table of clusters (idents)
  dat <-
    dplyr::left_join(clusters, meta.bc.group, by = "seurat_clusters")
  dat[is.na(dat)] = 0

  # convert to proportion
  dat <- dplyr::mutate(dat, Percentage = 100 * (.data$n / sum(.data$n)))

  # plot data
  if (plot) {
    if (isTRUE(plot.pct)) {
      p <-
        ggplot2::ggplot(dat, ggplot2::aes(
          x = !!idents,
          y = .data$Percentage,
          fill = !!idents
        )) +
        ggplot2::geom_histogram(stat = "identity") +
        ggplot2::scale_fill_manual(values = scales::hue_pal()(nrow(clusters))) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste("Percentage of cells in clusters:", group, factor))
    } else {
      p <- ggplot2::ggplot(dat, ggplot2::aes(
        x = !!idents,
        y = .data$n,
        fill = !!idents
      )) +
        ggplot2::geom_histogram(stat = "identity") +
        ggplot2::scale_fill_manual(values = scales::hue_pal()(nrow(clusters))) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste("Number of cells in clusters:", group, factor))
    }
    return(p)
  } else {
    return(dat)
  }
}
