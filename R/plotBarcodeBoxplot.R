#' @title
#' Boxplot of selected barcodes
#'
#' @description
#' Plots boxplot of counts for selected barcodes from DGEList object.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param barcodes Barcodes to be plotted (vector of strings).
#' @param group Optional, column name in sample metadata to facet data on (string).
#' @param conditions Optional, specific levels of group to be plotted (vector of strings).
#' @param trans Optional, transformation of y-axis, e.g. `log10` (string).
#' @param point Whether to include points (boolean). Default = `FALSE`.
#' @param violin Whether to include violin plots in addition to box plots (boolean). Default = `FALSE`.
#' @param returnData Whether to return data instead of plot (boolean). Default = `FALSE`.
#'
#' @return Returns a boxplot
#' @export
#' @examples
#' data(test.dge)
#' plotBarcodeBoxplot(test.dge, barcodes = "BC_190202", group = "Treatment",
#' conditions = c("Vehicle", "Low_dose", "High_dose"))

plotBarcodeBoxplot <- function(dgeObject,
                               barcodes = NULL,
                               group = NULL,
                               conditions = NULL,
                               trans = NULL,
                               point = FALSE,
                               violin = FALSE,
                               returnData = FALSE) {

  inputChecks(dgeObject, barcodes = barcodes, groups = group, conditions = conditions)

  counts.raw <- as.data.frame(dgeObject$counts)

  if (is.null(barcodes)) {
    stop("Please provide some barcodes to plot")
  }

  if (length(group) > 1) {
    stop("Please enter a single grouping variable")
  }

  if (!is.null(conditions)) {
    if (is.null(group)) {
      stop("Specify column to group by if selecting conditions.")
    }
  }

  sample = dgeObject$samples
  counts.cpm <- cpm(counts.raw)

  dat <-
    counts.cpm[which(rownames(counts.cpm) %in% barcodes), , drop = F]
  dat <- reshape2::melt(dat)
  colnames(dat) <- c("barcodes", "sample", "value")

  # get group data
  if (!is.null(group)) {
    group.dat <- dgeObject$samples[, group, drop = F]
    group.dat$sample <- rownames(group.dat)
    # incorporate metadata
    dat <- dplyr::left_join(dat, group.dat, by = "sample")
  }

  # subset on condition
  if (!is.null(conditions)) {
    dat <- dat[which(dat[, group] %in% conditions),]
  }
  if (returnData) {
    return(dat)
  }

  # plot logic:
  # if group and 1 barcode, then x = group
  # if group and >1 barcodes, then x=group, facets=barcodes
  # if no group and 1 barcode, then x=barcode
  # if no group and >1 barcode, then x=barcode

  # no metadata group specified
  if (is.null(group)) {
    p <- ggplot2::ggplot(dat) +
      ggplot2::geom_boxplot(aes(x = barcodes, y = value, fill = barcodes))

    if (violin) {
      p <- ggplot2::ggplot(dat) +
        ggplot2::geom_violin(aes(x = barcodes,
                                 y = value,
                                 fill = barcodes), scale = "width") +
        ggplot2::geom_boxplot(aes(x = barcodes, y = value), width = 0.2)
    }

    p <- p +
      ggplot2::xlab(NULL) +
      ggplot2::labs(xlab = as.factor(barcodes))

    if (point) {
      p <-
        p + ggplot2::geom_point(aes(x = barcodes, y = value))
    }

  } else {
    # check violin
    if (violin) {
      p <- ggplot2::ggplot(dat) +
        ggplot2::geom_violin(aes_string(y = "value",
                                        x = group,
                                        fill = group), scale = "width") +
        ggplot2::geom_boxplot(aes_string(y = "value", x = group), width = 0.2)
    } else {
      p <- ggplot2::ggplot(dat) +
        ggplot2::geom_boxplot(aes_string(y = "value",
                                         x = group,
                                         fill = group))
    }

    if (point) {
      p <- p +
        ggplot2::geom_point(aes_string(y = "value",
                                       x = group,
                                       fill = group))
    } # check point

    if (length(barcodes) > 1) {
      p <- p + ggplot2::labs(fill = as.character(group)) +
        facet_wrap(~ barcodes)
    }

  }
  p <- p +
    ggplot2::ylab("CPM") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(trans)) {
    p <- p + ggplot2::scale_y_continuous(trans = trans) +
      ggplot2::ylab(paste(trans, "CPM"))
  }

  return(p)
}
