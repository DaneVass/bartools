#'
#' @title
#' Plot barcode heatmap
#'
#' @description
#' Takes a DGEList with barcode counts, selects the n most abundant barcodes per sample and plots a heatmap.
#' Stars indicate the most abundant barcodes in the respective sample.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param topN Number of top barcodes per sample to show (integer). Default = `5`.
#' @param name String of the heatmap scale name (string). Default = `CPM`.
#' @param showBarcodes Show barcode names on the heatmap (boolean). Default = `FALSE`.
#' @param group Optional, column name in sample metadata to group samples by (string).
#' @param colAnnot List of vectors assigning colours to each level of metadata.
#' @param discrete Show presence/absence of barcodes instead of abundance (boolean). Default = `FALSE`.
#' @param discreteThreshold Threshold for presence of barcode in samples (decimal). Default = `1`.
#'
#' @return Returns a heatmap
#' @export
#' @examples
#' data(test.dge)
#' plotBarcodeHeatmap(test.dge)

plotBarcodeHeatmap <-
  function(dgeObject,
           topN = 5,
           name = "CPM",
           showBarcodes = FALSE,
           group = NULL,
           colAnnot = NULL,
           discrete = FALSE,
           discreteThreshold = 1) {
    if (methods::is(dgeObject)[1] != "DGEList") {
      stop("Please supply a valid DGEList object as input")
    }
    counts <- dgeObject$counts
    samples <- dgeObject$samples

    set.seed(1)
    bc <- c()
    # select top n barcodes for each sample
    for (s in colnames(counts)) {
      top_bc = utils::head(rownames(counts[order(counts[, s], decreasing = T),]), topN)
      bc <- c(bc, top_bc)
    }
    tab = counts[rownames(counts) %in% bc,]
    # add annotation of samples
    if (!is.null(group)) {
      # check that metadata names are in samples
      if (!all(group %in% colnames(samples))) {
        stop("group must be column in samples")
      }
      # keep group order for legend in plot
      group_cols <- group[group %in% colnames(samples)]
      dff <- samples[, group_cols, drop = FALSE]
      colnames(dff) <- group_cols

      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = dff,
        which = 'col',
        col = colAnnot,
        annotation_width = unit(c(1, 4), 'cm'),
        gap = unit(1, 'mm')
      )
    } else {
      ha <- NULL
    }
    if (discrete) {
      tab2 <-
        ifelse(counts[rownames(counts) %in% bc, ] > discreteThreshold, 1, 0)
      col <- c("1" = "red", "0" = "blue")

      suppressMessages(
        ComplexHeatmap::Heatmap(
          as.matrix(tab2),
          name = "Presence",
          show_row_names = showBarcodes,
          cell_fun = function(j, i, x, y, w, h, fill) {
            if (tab[i, j] %in% utils::head(sort(tab[, j], decreasing = TRUE), n = topN)) {
              grid::grid.text("*", x, y, vjust = 0.8)
            }
          },
          top_annotation = ha,
          col = col,
          heatmap_legend_param = list(labels = c("Yes", "No"))
        )
      )
    } else{
      suppressMessages(
        ComplexHeatmap::Heatmap(
          as.matrix(tab),
          name = name,
          show_row_names = showBarcodes,
          cell_fun = function(j, i, x, y, w, h, fill) {
            if (tab[i, j] %in% utils::head(sort(tab[, j], decreasing = TRUE), n = topN)) {
              grid::grid.text("*", x, y, vjust = 0.8)
            }
          },
          top_annotation = ha
        )
      )
    }
  }
