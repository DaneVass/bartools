#' plotBarcodesStackedBar
#'
#' Generate proportional stacked bar plot of barcodes from raw count object with n most frequent barcodes labelled
#'
#' @param counts dataframe containing raw counts of barcodes. Expects barcodes as row names and samples as columns. Alternatively, DGE object.
#' @param n_barcodes number of most frequent barcodes to show in legend
#' @param seed_colors seed for sampling colors
#' @param samples samples to plot
#'
#' @return Returns a stacked bar plot of barcode frequencies within samples
#'
#' @export
#'
#' @examples
#' plotBarcodesStackedBar(test.dge)

plotBarcodesStackedBar <-
  function(counts,
           n_barcodes = 10,
           seed_colors = 1,
           samples = NULL) {
    if (methods::is(counts)[1] == "DGEList") {
      # if no samplesheet is provided, extract it from DGE object
      counts <- as.data.frame(counts$counts)
    } else {
      counts <- as.data.frame(counts)
    }

    if (!is.null(samples)) {
      counts <- counts[, samples]
    }

    # calculate barcode frequencies within samples
    barcode_freqs <-
      counts %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "barcode") %>%
      pivot_longer(-barcode) %>%
      rename("count" = value, "sample" = name) %>%
      group_by(sample) %>%
      mutate(freq = count / sum(count)) %>%
      ungroup()

    # order barcodes based on maximum frequency across samples
    barcode_order <- barcode_freqs %>%
      select(barcode, freq) %>%
      slice_max(
        order_by = freq,
        by = barcode,
        n = 1,
        with_ties = F
      ) %>%
      arrange(freq)

    # set order of barcodes for plotting as factor levels
    barcode_freqs$barcode <- factor(barcode_freqs$barcode,
                                    levels = barcode_order$barcode)

    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
    col_vector = unlist(mapply(
      brewer.pal,
      qual_col_pals$maxcolors,
      rownames(qual_col_pals)
    ))

    n <- length(unique(barcode_freqs$barcode))

    # map colors to barcodes
    set.seed(seed_colors)
    cols_vector <- sample(col_vector, n, replace = T)

    # set which barcodes to label
    bc_labels <- barcode_freqs %>%
      arrange(desc(freq)) %>%
      pull(barcode) %>%
      unique() %>%
      head(n_barcodes) %>%
      as.character()

    label_values <-
      setNames(as.list(cols_vector), unique(barcode_freqs$barcode))

    p <- barcode_freqs %>%
      ggplot(aes(fill = barcode, x = sample, y = freq)) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = label_values, breaks = bc_labels) +
      theme_bw() +
      ylab("Proportion") +
      theme(legend.text = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1))

    return(p)
  }
