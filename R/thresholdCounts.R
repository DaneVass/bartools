#' Threshold counts
#'
#' Filter barcodes meeting a given absolute (total read count) or relative (proportion based) abundance level
#' Optionally plot number of barcodes detected using this threshold in each sample.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param type Threshold type to use (string). Must be one of "absolute" or "relative". Default = `absolute`.
#' @param threshold The threshold to use. If type = "relative", must a float between 0 & 1. Default = `10`
#' @param minSamples Minimum number of samples a barcode must meet threshold to remain in dataset (integer). Default = `1`.
#' @param plot Return plot instead of filtered data (boolean). Default = `FALSE`.
#' @param group Optional, column name in sample metadata to color samples by (string).
#' @param order Order samples by group (boolean). Default = `TRUE`.
#' @return Returns a filtered DGEList object.
#'
#' @export
#' @examples
#' data(test.dge)
#' thresholdCounts(test.dge, type = "absolute", threshold = 10)
#' thresholdCounts(test.dge, type = "absolute", threshold = 10, plot = TRUE)

thresholdCounts <- function(dgeObject,
                            threshold = 10,
                            type = "absolute",
                            minSamples = 1,
                            plot = FALSE,
                            group = NULL,
                            order = TRUE) {
  inputChecks(dgeObject, groups = group)
  # check obj

  # check threshold type
  if (!type %in% c("absolute", "relative")) {
    stop("type must be one of 'absolute' or 'relative'")
  }

  # check minSamples
  if (minSamples < 0 | minSamples > ncol(dgeObject)) {
    stop("minSamples must be greater than 0 and less or equal to ncol(dgeObject)")
  }

  # check threshold values
  if (type == "relative" &&
      threshold < 0 |
      type == "relative" &&
      threshold > 1 | !is.double(threshold)) {
    stop("relative threshold value must be a float between 0 and 1")
  }

  if (type == "absolute" && threshold < 0) {
    stop("absolute threshold value must be an integer > 0")
  }

  # check threshold for relative
  if (!is.null(threshold)) {
    threshold <- threshold
  } else {
    message("No threshold given, defaulting to type = 'absolute' & threshold = 10")
    type <- "absolute"
    threshold <- 10
  }

  # report dimensions pre and post
  message("DGEList dimensions pre-threshold")
  print(dim(dgeObject$counts))

  # filter DGEList based on thresholds
  if (type == "absolute") {
    keeprows = rowSums(dgeObject$counts >= threshold) >= as.numeric(minSamples)
    dgeObject <- dgeObject[keeprows, ]
  }

  if (type == "relative") {
    # convert everything to a proportion
    barcodes.proportional <- as.data.frame(dgeObject$counts)
    barcodes.proportional <-
      sweep(barcodes.proportional,
            2,
            colSums(barcodes.proportional),
            `/`)
    keeprows = rowSums(barcodes.proportional >= threshold) >= as.numeric(minSamples)
    dgeObject <- dgeObject[keeprows, ]
  }

  message("DGEList dimensions post-threshold")
  print(dim(dgeObject$counts))

  if (plot == FALSE) {
    # add number of detected barcodes above threshold to sample metadata
    above.threshold.counts <-
      data.frame(Sample = factor(), BC.count = c())

    for (sample in colnames(dgeObject$counts)) {
      above.threshold = length(which(dgeObject$counts[, sample] >= threshold))
      d <-
        data.frame(Sample = factor(sample), BC.count = above.threshold)
      above.threshold.counts <- rbind(above.threshold.counts, d)
    }

    # add metadata to object
    dgeObject$samples <- cbind(dgeObject$samples, above.threshold.counts)
    return(dgeObject)

  } else {
    above.threshold.counts <-
      data.frame(Sample = factor(), BC.count = c())

    for (sample in colnames(dgeObject$counts)) {
      above.threshold = length(which(dgeObject$counts[, sample] >= threshold))
      d <-
        data.frame(Sample = factor(sample), BC.count = above.threshold)
      above.threshold.counts <- rbind(above.threshold.counts, d)
    }

    if (is.null(group)) {
      g <- ggplot2::ggplot(above.threshold.counts,
                           ggplot2::aes(x = `Sample`, y = `BC.count`)) +
        ggplot2::geom_bar(stat = "identity")
    } else {
      above.threshold.counts$group <-
        dgeObject$samples[[group]]

      above.threshold.counts$group <-
        as.factor(above.threshold.counts$group)

      if (order) {
        above.threshold.counts$Sample <-
          factor(above.threshold.counts$Sample,
                 levels = above.threshold.counts$Sample[order(above.threshold.counts$group,
                                                              decreasing = T)])
      }

      g <- ggplot2::ggplot(above.threshold.counts,
                           ggplot2::aes(x = `Sample`, y = `BC.count`, fill = group)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(
          unique(above.threshold.counts$group)
        ))))
      g <- g + labs(fill = group)
    }
    g <- g +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour =
                                                                  "grey70")) +
      ggplot2::labs(
        title = paste(
          dim(dgeObject$counts)[[1]][1],
          " Total barcodes above ",
          type,
          " threshold = ",
          threshold,
          " in at least ",
          minSamples,
          " samples",
          sep = ""
        )
      ) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("Number of barcodes") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    return(g)
  }
}
