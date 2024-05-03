#' @title
#' Plot Barcode Similarity
#'
#' @description
#' Plot clustered barcode similarity heatmap based on detection across cells.
#'
#' @param bcSimilarity Similarity matrix with barcodes as rows and columns (output of `calculateBarcodeSimilarities()`).
#' @param topN Number of barcodes to display in heatmap, ranked by highest similarity metric (integer). Default = `50`.
#'
#' @return Returns a heatmap of barcode similarities.
#' @export
#'

plotBarcodeSimilarity <-
  function(bcSimilarity,
           topN = 50,
           dendrogram = F) {
    # get barcodes with highest similarity metric
    rowmax <- apply(bcSimilarity, 1, max, na.rm = TRUE)
    high_simil_bc <-
      names(head(rowmax[order(unlist(rowmax), decreasing = T)], n = topN))

    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(255)

    if (dendrogram) {
      # plot clustered heatmap
      p <- pheatmap::pheatmap(
        bcSimilarity[high_simil_bc, high_simil_bc],
        cluster_rows = T,
        cluster_cols = T,
        col = colors,
        display_numbers = F,
        cellheight = 15,
        cellwidth = 15,
        treeheight_row = 15,
        treeheight_col = 15,
        main = "",
      )
    } else {
      heatmap <-
        pheatmap::pheatmap(
          bcSimilarity[high_simil_bc, high_simil_bc],
          cluster_rows = T,
          cluster_cols = T,
          silent = T
        )
      roworder <- heatmap$tree_row$order

      # plot clustered heatmap
      p <- pheatmap::pheatmap(
        bcSimilarity[high_simil_bc, high_simil_bc][roworder, roworder],
        cluster_rows = F,
        cluster_cols = F,
        col = colors,
        display_numbers = F,
        cellheight = 15,
        cellwidth = 15,
        treeheight_row = 15,
        treeheight_col = 15,
        main = "",
      )
    }

    return(p)
  }
