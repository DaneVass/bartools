#' plotOrderedBubble
#'
#' Generate ordered proportional bubbleplots from raw count object with barcodes labelled above a specified threshold
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param orderSample Name of sample to order by (string).
#' @param proportionCutoff barcodes represented at a percentage within any sample above this threshold will be labelled (decimal). Default = `10`.
#' @param labelBarcodes Label barcodes with a proportion larger than proportionCutoff in any sample (boolean). Default = `TRUE`.
#' @param title Plot title (string). Default = `Proportional Bubble Plot`.
#' @param group Optional, column name in sample metadata to group samples by (string).
#' @param displaySamples Optional, vector of samples to display - keep the order of vector.
#' @param displayBarcodes Optional, vector of barcodes to display.
#' @param colorDominant Only color clones with frequency above `proportionCutoff` and others grey (boolean). Default = `FALSE`.
#' @param filterCutoff Barcodes below this threshold in `orderSample` will be filtered in all samples (boolean). Default = `TRUE`.
#' @param pseudoCount Whether to add a pseudo count of 1 to all counts to display barcodes absent in T0 (boolean). Requires counts to be normalized. Default = `FALSE`.
#' @param legend Show a legend of bubble sizes (boolean). Default = `TRUE`.
#'
#' @return Returns a bubbleplot of barcodes represented by proportion of total pool
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export
#' @examples
#' data(test.dge)
#' plotOrderedBubble(test.dge, orderSample = "T0-1", filterCutoff = 0.001, group = "Treatment")

plotOrderedBubble <- function(dgeObject,
                              title = "Proportional Bubble Plot",
                              orderSample = NULL,
                              group = NULL,
                              displaySamples = NULL,
                              displayBarcodes = NULL,
                              proportionCutoff = 10,
                              colorDominant = FALSE,
                              filterCutoff = NULL,
                              labelBarcodes = TRUE,
                              legend = TRUE,
                              pseudoCount = FALSE) {
  ###### check inputs ##########
  inputChecks(dgeObject, groups = group, samples = c(orderSample, displaySamples), barcodes = displayBarcodes)

  counts <- as.data.frame(dgeObject$counts)
  samples <- as.data.frame(dgeObject$samples)

  # check parameters

  if (is.null(orderSample)) {
    stop("Please provide sample to order by")
  }

  if (!is.null(group) & !all(group %in% colnames(samples))) {
    stop("Group must be column in samples")
  }

  # this will avoid plotting any barcode labels
  if (labelBarcodes == FALSE) {
    proportionCutoff = 100
  }

  if (pseudoCount == TRUE) {
    counts <- counts + 1
  }

  # transform CPM into percentage within sample
  barcodes.proportional <-
    sweep(counts,
          2,
          colSums(counts),
          `/`) * 100
  if (!is.null(displaySamples)) {
    # filter and order based on displaySamples
    barcodes.proportional <-
      barcodes.proportional[, displaySamples]
  }

  # Order by selected sample.
  barcodes.proportional$Position <-
    barcodes.proportional[, orderSample]
  barcodes.proportional$Barcode <- rownames(barcodes.proportional)
  if (!is.null(filterCutoff)) {
    barcodes.proportional <-
      barcodes.proportional[barcodes.proportional$Position > filterCutoff,]
  }
  # identify all barcodes to label
  Highbarcodes <- barcodes.proportional %>%
    dplyr::select(-c(.data$Barcode, .data$Position)) %>%
    dplyr::filter(dplyr::if_any(tidyselect::where(is.numeric), ~ .x > proportionCutoff))

  if (colorDominant) {
    # make all barcodes grey and only color those that are above threshold cutoff
    colors <- "grey90"
  } else {
    # give each barcode a specific color
    colors <- scales::hue_pal()(30)
  }

  # filter barcodes to display before sampling colors
  if (!is.null(displayBarcodes)) {
    barcodes.proportional <-
      barcodes.proportional[rownames(barcodes.proportional) %in% displayBarcodes, ]
    Highbarcodes <-
      Highbarcodes[rownames(Highbarcodes) %in% displayBarcodes,]
  }

  colors <-
    sample(colors, length(rownames(barcodes.proportional)), replace = TRUE)
  names(colors) <- rownames(barcodes.proportional)
  barcodes.proportional$Color <- colors

  if (nrow(Highbarcodes) != 0) {
    # Assign diverse colours to high frequency barcodes
    SelColors <- scales::hue_pal()(nrow(Highbarcodes))
    i = 1
    for (bc in rownames(Highbarcodes)) {
      barcodes.proportional[bc,]$Color <- SelColors[i]
      i <- i + 1
    }
  }
  HighbarcodesLabel <-
    barcodes.proportional[rownames(Highbarcodes),]
  HighbarcodesLabel <-
    HighbarcodesLabel[HighbarcodesLabel$Position > 0,]

  # melt data frame and rename columns correctly
  barcodes.proportional.melted <-
    reshape2::melt(barcodes.proportional,
                   id.vars = c("Color", "Position", "Barcode"))
  colnames(barcodes.proportional.melted) <-
    c("Color", "Position", "Barcode", "Sample", "Proportion")
  # convert variables to correct form
  barcodes.proportional.melted$Barcode <-
    as.factor(barcodes.proportional.melted$Barcode)
  # barcodes.proportional.melted$Sample <-
  #   as.factor(barcodes.proportional.melted$Sample)
  # ordered sample on top by ordering levels of factor
  barcodes.proportional.melted$Sample <-
    factor(
      barcodes.proportional.melted$Sample,
      levels = c(as.character(
        unique(barcodes.proportional.melted$Sample[barcodes.proportional.melted$Sample != orderSample])
      ), orderSample),
      ordered = TRUE
    )
  barcodes.proportional.melted$Proportion <-
    as.numeric(barcodes.proportional.melted$Proportion)

  # add metadata to bubbleplot
  if (!is.null(group)) {
    # check if all samples should be displayed
    if (!is.null(displaySamples)) {
      samples.filtered <- data.frame(samples[displaySamples, group])
      samples.filtered$Sample <- displaySamples
    } else {
      samples.filtered <- data.frame(samples[, group])
      samples.filtered$Sample <- colnames(counts)
    }
    colnames(samples.filtered) <- c(group, "Sample")

    # order groups so that group containing ordered sample is on top
    orderGroup <-
      samples.filtered[samples.filtered$Sample == orderSample, group]
    samples.filtered[group] <-
      factor(
        samples.filtered[[group]],
        levels = c(orderGroup, as.character(unique(
          samples.filtered[group][samples.filtered[group] != orderGroup]
        ))),
        ordered = TRUE
      )

    # create new merged data frame which including barcode data with metadata for the samples
    barcodes.proportional.melted <-
      merge(
        barcodes.proportional.melted,
        samples.filtered,
        by = "Sample",
        all.x = TRUE
      )
  }

  # generate bubbleplot
  bubble.plot <- ggplot2::ggplot(
    barcodes.proportional.melted,
    ggplot2::aes(
      x = .data$Position,
      y = .data$Sample,
      size = .data$Proportion,
      color = .data$Color
    )
  ) +
    ggplot2::geom_point(stat = "identity",
                        alpha = 0.6,
                        shape = 16) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(y = "Condition",
                  x = "",
                  title = title) +
    ggplot2::scale_size_continuous(
      name = "Barcode\nProportion (%)",
      range = c(0.1, 10),
      breaks = c(0.1, 1, 2, 5, 10, 20, 40, 60, 80),
      labels = c(0.1, 1, 2, 5, 10, 20, 40, 60, 80),
    ) +
    ggplot2::scale_x_continuous(
      trans = 'log10',
      labels = HighbarcodesLabel$Barcode,
      breaks = HighbarcodesLabel$Position,
      sec.axis = ggplot2::sec_axis(
        ~ . * 1,
        labels = c(0.0001, 0.001, 0.01, 0.1, 1, 2, 5, 10, 20, 30, 40, 50),
        breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 2, 5, 10, 20, 30, 40, 50),
        name = paste("Barcode Proportion in", orderSample, "(%)")
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 6),
      legend.box.spacing = unit(2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0, "mm"),
      legend.spacing.y = unit(0, "mm"),
      axis.text.y = ggplot2::element_text(size = 6),
      plot.title = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 6)
    )
  # turn and color labels if labels are plotted
  if (nrow(Highbarcodes) > 0) {
    bubble.plot <- bubble.plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0.5,
          colour = HighbarcodesLabel$Color,
          size = 6
        )
      )
  }

  # facet plot if grouping is provided
  if (!is.null(group)) {
    bubble.plot <- bubble.plot +
      facet_grid(barcodes.proportional.melted[[group]] ~ .,
                 scales = "free",
                 space = "free")
  }

  # remove legend
  if (legend == FALSE) {
    bubble.plot <- bubble.plot +
      theme(legend.position = "None")
  }

  return(bubble.plot)
}
