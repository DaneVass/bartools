#' plotDetectedBarcodes
#'
#' Plot the number of barcodes comprising the top Nth percentile for each sample as a bar plot.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param percentile Percentile threshold to count barcodes (decimal). Default = `0.95` (i.e. 95th percentile).
#' @param plot Plot data instead of returning counts table (boolean).
#' Table of number of barcodes and lists of barcodes in the top Nth percentile can be obtained with `calcPercentileBarcodes()`. Default = `TRUE`.
#' @param sampleOrder Optional, Ordering of samples (vector of strings).
#' @param group Optional, column name in sample metadata to group samples by (string).
#' @param title Optional, Plot title (string).
#'
#' @return Returns a histogram plot of the number of detected barcodes per sample in the top Nth percentile.
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' data(test.dge)
#' plotDetectedBarcodes(test.dge, percentile = .95)

plotDetectedBarcodes <-
  function(dgeObject,
           percentile = 0.95,
           plot = T,
           sampleOrder = NULL,
           group = NULL,
           title = NULL) {
    ###### check inputs ##########
    inputChecks(dgeObject, groups = group, samples = sampleOrder)

    # reorder or subset DGEList object if sample order is given
    if (is.null(sampleOrder)) {
      sampleOrder <- colnames(dgeObject)
    } else {
      dgeObject <- dgeObject[, sampleOrder]
    }

    if (!is.null(group)) {
      cols <- as.factor(group)
    }
    else {
      cols <- rep("Group 1", length(colnames(dgeObject)))
    }

    # calculate number of barcodes in the top Nth percentile per sample
    percentile.df <-
      calcPercentileBarcodes(dgeObject, percentile = percentile)$NumBarcodes

    # merge group into percentile.df
    if (!is.null(group)) {
      groups <-
        dgeObject$samples[, which(colnames(dgeObject$samples) == as.character(group)), drop = F]
      groups$Sample <- rownames(groups)
      groups <- groups[, c(2, 1)]
      colnames(groups) <- c("Sample", "Group")
      percentile.df <-
        suppressMessages(dplyr::left_join(percentile.df, groups))
      percentile.df$Sample <-
        factor(percentile.df$Sample, levels = sampleOrder)
      percentile.df$Group <- as.factor(percentile.df$Group)
      percentile.df <-
        dplyr::mutate(percentile.df,
                      Sample = forcats::fct_reorder(.data$Sample, dplyr::desc(.data$Group)))
    }

    # setup plot title
    if (is.null(title)) {
      title <-
        paste("Number of detected barcodes per sample. Percentile =",
              as.character(percentile))
    }

    # make plots, color by group if given, needs ggpubr
    if (plot) {
      if (!is.null(group)) {
        p <-
          ggplot2::ggplot(data = percentile.df,
                          ggplot2::aes(
                            x = .data$Sample,
                            y = .data$NumBarcodes,
                            fill = .data$Group
                          ))
      } else {
        p <-
          ggplot2::ggplot(data = percentile.df,
                          ggplot2::aes(x = .data$Sample,
                                       y = .data$NumBarcodes))
      }
      p <- p + ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0,
          hjust = 1
        )) + ggplot2::ggtitle(title) +
        ggplot2::ylim(0, max(percentile.df$NumBarcodes + 25)) +
        ggplot2::geom_text(
          ggplot2::aes(
            label = .data$NumBarcodes,
            y = .data$NumBarcodes + 0.05
          ),
          position = ggplot2::position_dodge(0.9),
          vjust = -0.1
        ) +
        ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(
          unique(percentile.df$Group)
        ))))
      return(p)
    }
    else {
      return(percentile.df)
    }
  }
