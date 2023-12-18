#' Plot sample read counts
#'
#' Simple plot of total read counts per sample
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param group Optional, column name in sample metadata to group samples by (string).
#' @param log10 log10 transform output (boolean). Default = `FALSE`.
#' @param legend Include legend (boolean). Default = `TRUE`.
#' @param order Order samples by group (boolean). Default = `TRUE`.
#'
#' @return Returns a plot of the read counts per column (sample) in a data frame
#' @export
#'
#' @examples
#' data(test.dge)
#' plotReadCounts(test.dge, group = "Treatment")
#'

plotReadCounts <-
  function(dgeObject,
           group = NULL,
           log10 = FALSE,
           legend = TRUE,
           order = TRUE) {
    inputChecks(dgeObject, groups = group)

    counts <- dgeObject$counts

    if (!is.null(group)) {
      cols <- as.factor(dgeObject$samples[[group]])
    } else {
      cols <- rep("Group 1", length(colnames(counts)))
    }

    dat <- colSums(counts)
    dat <-
      data.frame(
        sample = names(dat),
        counts = dat,
        group = cols,
        row.names = NULL
      )

    if (order) {
      dat$sample <-
        factor(dat$sample, levels = dat$sample[order(dat$group, decreasing = T)])
    }

    # plot data
    p <-
      ggplot2::ggplot(data = dat, ggplot2::aes(y = sample, x = counts, fill = group)) +
      ggplot2::geom_bar(stat = "identity", width = .75) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, vjust =
                                                           0.5)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5)) +
      ggplot2::ggtitle("Total Read Counts") +
      ggplot2::xlim(0, max(dat$counts + 100)) +
      ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(
        unique(dat$group)
      ))))

    if (isTRUE(log10)) {
      p <- p + ggplot2::scale_x_log10()
    }

    if (!is.null(group)) {
      p <- p + labs(fill = group)
    } else {
      legend = FALSE
    }

    if (legend == FALSE) {
      p <- p + ggplot2::theme(legend.position = "none")
    }

    return(p)
  }
