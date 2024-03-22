#' @title
#' Calculate barcode similarities
#'
#' @description
#' Calculate barcode similarities in single-cell data in order to group barcodes into clones.
#'
#' @param counts Count table with columns cellid, barcode and bc.umi.count (output of `readBartabCounts()`).
#' @param umiCutoff Minimum UMI count to consider a barcode detected in a cell (integer). Default = `1`.
#' @param metric Similarity metric to use, one of `jaccard` and `pearson` (string). Default = `jaccard`.
#' @param plotHistogram Whether to plot a histogram of the similarity scores (bool). Default = `TRUE`.
#'
#' @return Returns a similarity matrix of all barcodes.
#' @export
#'


calculateBarcodeSimilarities <-
  function(counts,
           umiCutoff = 1,
           metric = "jaccard",
           plotHistogram = TRUE) {
    if (!metric %in% c("jaccard", "pearson")) {
      stop("Metric must be one of jaccard, pearson")
    }
    # turn counts from readBartabCounts() into wide matrix
    counts_matrix <- counts %>%
      tidyr::pivot_wider(id_cols = cellid,
                         values_from = bc.umi.count,
                         names_from = barcode)

    counts_matrix <- as.data.frame(counts_matrix)
    rownames(counts_matrix) <- counts_matrix$cellid
    counts_matrix$cellid <- NULL

    # binarize count matrix given a UMI threshold
    counts_matrix_binary <- counts_matrix
    counts_matrix_binary[is.na(counts_matrix_binary)] <- 0
    counts_matrix_binary[counts_matrix_binary < umiCutoff] <- 0
    counts_matrix_binary[counts_matrix_binary >= umiCutoff] <- 1

    # calcualte pearson correlation
    if (metric == "pearson") {
      bc_similarity <- stats::cor(counts_matrix_binary)
      # set diagonal to NA
      diag(bc_similarity) <- NA
      hist_title <- "Histogram of pearson correlation coefficient"
    }
    # or calculate jaccard index
    else if (metric == "jaccard") {
      # proxyC is slightly rounded but definitely negligible effect.
      bc_similarity <-
        proxyC::simil(t(as.matrix(counts_matrix_binary)), method = "jaccard")
      bc_similarity <- as.matrix(bc_similarity)
      # set diagonal to NA
      diag(bc_similarity) <- NA
      hist_title <- "Histogram of jaccard index"
    }

    # plot histogram of similarity metric
    if (plotHistogram) {
      hist(bc_similarity, breaks = 100, main = hist_title)
    }
    return(bc_similarity)
  }
