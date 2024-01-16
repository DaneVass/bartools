#' plotDetectedBarcodes
#'
#' Plot the total number of barcodes detected in a sample above a percentile threshold.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param percentile Percentile threshold to count barcodes (decimal). Default = `0.95`.
#' @param plot Plot data instead of returning counts table (boolean). Default = `TRUE`.
#' @param sampleOrder Optional, Ordering of samples (vector of strings).
#' @param group Optional, column name in sample metadata to group samples by (string).
#' @param title Optional, Plot title (string).
#'
#' @return Returns a histogram plot of the number of detected barcodes per sample.
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' data(test.dge)
#' plotDetectedBarcodes(test.dge, percentile = .95)
#' plotDetectedBarcodes(test.dge, plot = FALSE)

plotDetectedBarcodes <-
  function(dgeObject,
           percentile = 0.95,
           plot = T,
           sampleOrder = NULL,
           group = NULL,
           title = NULL) {
    ###### check inputs ##########
    inputChecks(dgeObject, groups = group, samples = sampleOrder)

    counts.obj <- as.data.frame(dgeObject$counts)

    if (!is.null(group)) {
      cols <- as.factor(group)
    }
    else {
      cols <- rep("Group 1", length(colnames(dgeObject)))
    }

    # reorder counts object if given
    if (is.null(sampleOrder)) {
      sampleOrder <- colnames(dgeObject)
    }
    counts.obj <- as.data.frame(counts.obj[, c(sampleOrder)])

    samples <- colnames(counts.obj)
    percentile.df <- data.frame()
    barcodes <- c()
    for (i in samples) {
      if (i == "color") {
        message("skipping color column")
      }
      else {
        order <- order(counts.obj[, i], decreasing = T)
        sorted <- counts.obj[order, as.character(i), drop = F]
        sorted <- sorted[sorted > 0, , drop = F]
        colsum <- sum(sorted)
        percentile.cutoff <- sum(sorted) * percentile
        cumsum <- cumsum(sorted)
        stopifnot(max(cumsum) == colsum)
        len <- length(which(cumsum <= percentile.cutoff))
        d <- data.frame(Sample = factor(i), Barcodes = len)
        percentile.df <- rbind(percentile.df, d)
        rows <- rownames(sorted)[which(cumsum <= percentile.cutoff)]
        length(rows)
        sorted.top <- counts.obj[rows, i, drop = F]
        top.bc.sum <- sum(sorted.top[, 1])
        sorted.top$percentage <-
          (sorted.top[, 1] / top.bc.sum) * 100
      }
    }

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
                            y = .data$Barcodes,
                            fill = .data$Group
                          ))
      } else {
        p <-
          ggplot2::ggplot(data = percentile.df, ggplot2::aes(x = .data$Sample,
                                                             y = .data$Barcodes))
      }
      p <- p + ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0,
          hjust = 1
        )) + ggplot2::ggtitle(title) +
        ggplot2::ylim(0, max(percentile.df$Barcodes + 25)) +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$Barcodes, y = .data$Barcodes + 0.05),
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
