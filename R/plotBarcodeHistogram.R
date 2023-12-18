#' plotBarcodeHistogram
#'
#' Generate proportional stacked bar plot of barcodes from raw count object with n most frequent barcodes labelled.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param topN Number of most frequent barcodes to show in legend (integer). Default = `10`.
#' @param seedColors Seed for sampling colors (integer). Default = `1`.
#' @param samples Samples to plot (vector of strings). Default all.
#' @param orderSamples One or multiple samples to order barcodes by (string or vector of strings). Default all.
#' @param alphaLowFreq Alpha of barcodes not in top n barcodes (decimal). Default = `1`.
#'
#' @return Returns a stacked bar plot of barcode frequencies within samples
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' data(test.dge)
#' plotBarcodeHistogram(test.dge)
plotBarcodeHistogram <-
  function(dgeObject,
           topN = 10,
           seedColors = 1,
           samples = NULL,
           orderSamples = NULL,
           alphaLowFreq = 1) {

    inputChecks(dgeObject, samples = c(samples, orderSamples))

    if (topN > 74) {
      message("Warning: number of labelled barcodes larger than number of colors in palette (74)")
    }

    counts <- as.data.frame(dgeObject$counts)

    if (is.null(samples)) {
      samples <- colnames(counts)
    }

    counts <- counts[, samples]

    if (!is.null(orderSamples)) {
      if (!all(orderSamples %in% samples)) {
        stop("samples to order by are not present or among selected samples")
      }
    } else {
      orderSamples <- samples
    }

    # calculate barcode frequencies within samples
    barcode_freqs <-
      counts %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "barcode") %>%
      tidyr::pivot_longer(-.data$barcode) %>%
      dplyr::rename("count" = .data$value, "sample" = .data$name) %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::mutate(freq = .data$count / sum(.data$count)) %>%
      dplyr::ungroup()

    # order barcodes based on maximum frequency across samples
    barcode_order <- barcode_freqs %>%
      dplyr::filter(.data$sample %in% orderSamples) %>%
      dplyr::select(.data$barcode, .data$freq) %>%
      dplyr::group_by(.data$barcode) %>%
      dplyr::slice_max(
        order_by = .data$freq,
        n = 1,
        with_ties = F
      ) %>%
      dplyr::arrange(.data$freq)

    # set order of barcodes for plotting as factor levels
    barcode_freqs$barcode <- factor(barcode_freqs$barcode,
                                    levels = barcode_order$barcode)

    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
    col_vector = unlist(mapply(
      RColorBrewer::brewer.pal,
      qual_col_pals$maxcolors,
      rownames(qual_col_pals)
    ))

    # TODO make sure that topN barcodes have different colors -> sample separately
    n <- length(unique(barcode_freqs$barcode))

    # map colors to barcodes
    set.seed(seedColors)
    cols_vector <- sample(col_vector, n, replace = T)

    # set which barcodes to label
    bc_labels <- barcode_freqs %>%
      dplyr::arrange(dplyr::desc(.data$freq)) %>%
      dplyr::pull(.data$barcode) %>%
      unique() %>%
      utils::head(topN) %>%
      as.character()

    label_values <-
      stats::setNames(as.list(cols_vector), unique(barcode_freqs$barcode))

    p <- barcode_freqs %>%
      ggplot2::ggplot(aes(
        fill = .data$barcode,
        x = .data$sample,
        y = .data$freq,
        alpha = ifelse(.data$barcode %in% bc_labels, 1, 0.5)
      )) +
      ggplot2::geom_bar(position = "stack", stat = "identity") +
      ggplot2::scale_alpha_continuous(guide = FALSE, range = c(alphaLowFreq, 1)) +
      ggplot2::scale_fill_manual(values = label_values, breaks = bc_labels) +
      ggplot2::theme_bw() +
      ggplot2::labs(fill = paste0("Top ", topN, " barcodes")) +
      ggplot2::ylab("Proportion") +
      ggplot2::theme(legend.text = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1))

    return(p)
  }
