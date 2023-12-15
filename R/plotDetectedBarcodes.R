#' plotDetectedBarcodes
#'
#' Plot the total number of barcodes detected in a sample
#'
#' @param counts DGEList or dataframe containing raw or normalised barcode counts
#' @param percentile desired percentile value. 95th percentile by default
#' @param plot Logical. plot data instead of returning counts table.
#' @param sample.order desired ordering of the samples on the plot
#' @param group grouping field in dgelist$samples to color samples by
#' @param title desired plot title
#'
#' @return Returns a histogram plot of the number of detected barcodes per sample.
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' plotDetectedBarcodes(test.dge, percentile = .95)
#' plotDetectedBarcodes(test.dge, plot = FALSE)

plotDetectedBarcodes <- function(counts, percentile = 0.95, plot = T, sample.order = NULL, group = NULL, title = NULL)
{
  # check inputs
  if (methods::is(counts)[1] == "DGEList") {
    counts.obj <- as.data.frame(counts$counts)
  }
  else {
    counts.obj <- as.data.frame(counts)
  }

  if (!is.null(group)) {
    cols <- as.factor(group)
  }
  else {
    cols <- rep("Group 1", length(colnames(counts)))
  }

  # reorder counts object if given
  if (is.null(sample.order)){
    counts.obj <- counts.obj
    #message("default sample ordering as per colnames(dge)")
  } else {
    if (length(sample.order) != length(colnames(counts.obj))){
      stop("sample.order must be a vector of same length as colnames(dge)")
    }
    if (sort(as.numeric(sample.order)) != 1:length(colnames(counts.obj))){
      stop("sample.order must contain the same values as colnames(dge)")
    }
    counts.obj <- as.data.frame(counts.obj[,c(sample.order)])
  }


  #dim(counts.obj)
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
      sorted.top$percentage <- (sorted.top[, 1]/top.bc.sum) * 100
    }
  }

  # merge group into percentile.df
  if (methods::is(counts)[1] == "DGEList") {
    if (!is.null(group)) {
      groups <- counts$samples[,which(colnames(counts$samples) == as.character(group)), drop = F]
      groups$Sample <- rownames(groups)
      groups <- groups[,c(2,1)]
      colnames(groups) <- c("Sample", "Group")
      percentile.df <- suppressMessages(dplyr::left_join(percentile.df, groups))
      percentile.df$Sample <- as.factor(percentile.df$Sample)
      percentile.df$Group <- as.factor(percentile.df$Group)
      percentile.df <- dplyr::mutate(percentile.df, Sample = forcats::fct_reorder(Sample, dplyr::desc(Group)))
    }
  }

  # setup plot title
  if (is.null(title)){
    title <- paste("Number of detected barcodes per sample. Percentile =", as.character(percentile))
  }

  # make plots, color by group if given, needs ggpubr
  if (plot) {
    if (!is.null(group)){
      p <- ggplot2::ggplot(data = percentile.df, ggplot2::aes(x = `Sample`, y = `Barcodes`, fill = Group)) +
        ggplot2::geom_bar(stat = "identity") +
        #ggplot2::guides(fill = "none") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0, hjust = 1, size = 5)) + ggplot2::ggtitle(title) +
        ggplot2::ylim(0, max(percentile.df$Barcodes + 25)) +
        ggplot2::geom_text(ggplot2::aes(label = `Barcodes`, y = `Barcodes` + 0.05), position = ggplot2::position_dodge(0.9), vjust = -0.1, size = 5) +
        ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(unique(percentile.df$Group)))))
      return(p)
    } else {
      p <- ggplot2::ggplot(data = percentile.df, ggplot2::aes(x = percentile.df$Sample,
                                                              y = percentile.df$Barcodes)) +
        ggplot2::geom_bar(stat = "identity") +
        #ggplot2::guides(fill = "none") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0, hjust = 1, size = 5)) + ggplot2::ggtitle(title) +
        ggplot2::ylim(0, max(percentile.df$Barcodes + 25)) +
        ggplot2::geom_text(ggplot2::aes(label = percentile.df$Barcodes, y = percentile.df$Barcodes + 0.05), position = ggplot2::position_dodge(0.9), vjust = -0.1, size = 5) +
        ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(unique(percentile.df$Group)))))
      return(p)
    }
  }
  else {
    return(percentile.df)
  }
}
