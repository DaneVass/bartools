#' @title
#' Boxplot of selected barcodes
#'
#' @description
#' Plots boxplot of counts for selected barcodes from DGEList object.
#'
#' @param dgeObject DGEList object with normalized barcode counts.
#' @param barcodes Barcodes to be plotted (vector of strings).
#' @param group Optional, column name in sample metadata to facet data on (string).
#' @param conditions Optional, specific levels of group to be plotted (vector of strings).
#' @param trans Optional, transformation of y-axis, e.g. `log10` (string).
#' @param point Whether to include points (boolean). Default = `FALSE`.
#' @param violin Whether to include violin plots in addition to box plots (boolean). Default = `FALSE`.
#' @param returnData Whether to return data instead of plot (boolean). Default = `FALSE`.
#' @param normalizeMethod Method for normalizing counts (string). One of `CPM`, `TMM`, `TMMwsp`, `upperquartile`, `RLE` or `NULL` for no normalization. See `edgeR::calcNormFactors()`. Default = `CPM`.
#'
#'
#' @return Returns a boxplot
#' @importFrom rlang .data
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
                               returnData = FALSE,
                               normalizeMethod = "CPM") {

  inputChecks(dgeObject, barcodes = barcodes, groups = group, conditions = conditions)

  if (!is.null(normalizeMethod) & normalizeMethod != "NULL") {
    ylabel <- normalizeMethod
    dgeObject <- normaliseCounts(dgeObject = dgeObject, method = normalizeMethod)
  } else {
    ylabel <- "counts"
  }

  counts <- dgeObject$counts

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

  dat <-
    counts[which(rownames(counts) %in% barcodes), , drop = F]
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
      ggplot2::geom_boxplot(aes(x = .data$barcodes, y = .data$value, fill = .data$barcodes))

    if (violin) {
      p <- ggplot2::ggplot(dat) +
        ggplot2::geom_violin(aes(x = .data$barcodes,
                                 y = .data$value,
                                 fill = .data$barcodes), scale = "width") +
        ggplot2::geom_boxplot(aes(x = .data$barcodes, y = .data$value), width = 0.2)
    }

    p <- p +
      ggplot2::xlab(NULL) +
      ggplot2::labs(xlab = as.factor(barcodes))

    if (point) {
      p <-
        p + ggplot2::geom_point(aes(x = .data$barcodes, y = .data$value))
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
        ggplot2::facet_wrap(~ .data$barcodes)
    }

  }
  p <- p +
    ggplot2::ylab(ylabel) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(trans)) {
    p <- p + ggplot2::scale_y_continuous(trans = trans) +
      ggplot2::ylab(paste(trans, ylabel))
  }

  return(p)
}
