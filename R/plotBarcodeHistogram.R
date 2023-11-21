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
#'
#' @export
#'
#' @examples
#' plotBarcodesStackedBar(test.dge)

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
      tidyr::pivot_longer(-barcode) %>%
      dplyr::rename("count" = value, "sample" = name) %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(freq = count / sum(count)) %>%
      dplyr::ungroup()

    # order barcodes based on maximum frequency across samples
    barcode_order <- barcode_freqs %>%
      dplyr::filter(sample %in% orderSamples) %>%
      dplyr::select(barcode, freq) %>%
      dplyr::group_by(barcode) %>%
      dplyr::slice_max(
        order_by = freq,
        n = 1,
        with_ties = F
      ) %>%
      dplyr::arrange(freq)

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
      dplyr::arrange(desc(freq)) %>%
      dplyr::pull(barcode) %>%
      unique() %>%
      head(topN) %>%
      as.character()

    label_values <-
      setNames(as.list(cols_vector), unique(barcode_freqs$barcode))

    p <- barcode_freqs %>%
      ggplot(aes(
        fill = barcode,
        x = sample,
        y = freq,
        alpha = ifelse(barcode %in% bc_labels, 1, 0.5)
      )) +
      geom_bar(position = "stack", stat = "identity") +
      scale_alpha_continuous(guide = FALSE, range = c(alphaLowFreq, 1)) +
      scale_fill_manual(values = label_values, breaks = bc_labels) +
      theme_bw() +
      labs(fill = paste0("Top ", topN, " barcodes")) +
      ylab("Proportion") +
      theme(legend.text = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1))

    return(p)
  }
