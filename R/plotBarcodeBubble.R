#' plotBarcodeBubble
#'
#' Generate proportional bubbleplots from raw count object with barcodes labelled above a specified threshold
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param proportionCutoff Barcodes represented at a percentage within any sample above this threshold will be labelled (decimal). Default = `10`.
#' @param labelBarcodes Label barcodes with a proportion larger than `proportionCutoff` in any sample (boolean). Default = `TRUE`.
#' @param title Plot title (string). Default = `Proportional Bubble Plot`.
#' @param group Optional, column name in sample metadata to group samples by (string).
#' @param displaySamples Optional, vector of sample names to display, preserves the order of vector (vector of strings).
#' @param displayBarcodes Optional, vector of barcodes to display (vector of strings).
#' @param colorDominant Only color clones with frequency above `proportionCutoff` and others colored grey (boolean). Default = `FALSE`.
#' @param legend Show a legend of bubble sizes (boolean). Default = `TRUE`.
#' @param orderBarcodes Order barcodes alphanumerical (boolean). For SPLINTR that represents abundance in  original barcode library. Default `TRUE`.
#'
#' @return Returns a bubbleplot of barcodes represented by proportion of total pool
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' data(test.dge)
#' plotBarcodeBubble(test.dge, proportionCutoff = 10)

plotBarcodeBubble <- function(dgeObject,
                              title = "Proportional Bubble Plot",
                              group = NULL,
                              displaySamples = NULL,
                              displayBarcodes = NULL,
                              proportionCutoff = 10,
                              colorDominant = FALSE,
                              labelBarcodes = TRUE,
                              orderBarcodes = TRUE,
                              legend = TRUE) {
  ###### check inputs ##########
  inputChecks(dgeObject, groups = group, samples = displaySamples, barcodes = displayBarcodes)

  counts <- dgeObject$counts
  samples <- dgeObject$samples

  # this will avoid plotting any barcode labels
  if (labelBarcodes == FALSE) {
    proportionCutoff = 100
  }

  # transform counts / CPM into percentage within sample
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

  # maintain Rank information of barcodes. Plot in ascending rank order from original barcode library.
  # this is done by ordering barcodes alphanumerical.
  # this only makes sense for SPLINTR, where the rank order is in the barcode name
  # but not barcode systems like DRAG and is thus optional but default behavior.
  if (orderBarcodes) {
    barcodes.proportional <-
      barcodes.proportional[stringr::str_sort(rownames(barcodes.proportional), numeric = T),]
    barcodes.proportional <- as.data.frame(barcodes.proportional)
  }
  barcodes.proportional$Position <-
    as.factor(seq(1, length(rownames(counts))))
  barcodes.proportional$Barcode <- rownames(barcodes.proportional)

  # identify all barcodes to label
  Highbarcodes <- barcodes.proportional %>%
    dplyr::select(-c(Barcode, Position)) %>%
    dplyr::filter(dplyr::if_any(tidyselect::where(is.numeric), ~ .x > proportionCutoff))

  if (colorDominant) {
    # make all barcodes grey and only color those that are above threshold cutoff
    colors <- "grey90"
  } else {
    # give each barcode a specific color
    colors <- scales::hue_pal()(30)
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
    HighbarcodesLabel[as.numeric(HighbarcodesLabel$Position) > 0,]

  if (!is.null(displayBarcodes)) {
    barcodes.proportional <-
      barcodes.proportional[rownames(barcodes.proportional) %in% displayBarcodes, ]
    HighbarcodesLabel <-
      HighbarcodesLabel[HighbarcodesLabel$Barcode %in% displayBarcodes,]
  }

  # melt data frame and rename columns correctly
  barcodes.proportional.melted <-
    reshape2::melt(barcodes.proportional,
                   id.vars = c("Color", "Position", "Barcode"))
  colnames(barcodes.proportional.melted) <-
    c("Color", "Position", "Barcode", "Sample", "Proportion")

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

    # create new merged data frame which including barcode data with metadata for the samples
    barcodes.proportional.melted <-
      merge(
        barcodes.proportional.melted,
        samples.filtered,
        by = "Sample",
        all.x = TRUE
      )
  }

  # convert variables to correct form
  barcodes.proportional.melted$Barcode <-
    as.factor(barcodes.proportional.melted$Barcode)
  barcodes.proportional.melted$Sample <-
    as.factor(barcodes.proportional.melted$Sample)
  barcodes.proportional.melted$Proportion <-
    as.numeric(barcodes.proportional.melted$Proportion)

  # generate bubbleplot
  bubble.plot <- ggplot2::ggplot(
    barcodes.proportional.melted,
    ggplot2::aes(
      x = Position,
      y = Sample,
      size = Proportion,
      color = Color
    )
  ) +
    ggplot2::geom_point(stat = "identity",
                        alpha = 0.6,
                        shape = 16) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(y = "Condition",
                  x = "Barcode",
                  title = title) +
    ggplot2::scale_size_continuous(
      name = "Barcode\nProportion (%)",
      range = c(0.1, 10),
      breaks = c(0.1, 1, 2, 5, 10, 20, 40, 60, 80),
      labels = c(0.1, 1, 2, 5, 10, 20, 40, 60, 80),
    ) +
    ggplot2::scale_x_discrete(breaks = HighbarcodesLabel$Position,
                              labels = HighbarcodesLabel$Barcode) +
    ggplot2::theme_bw()

  # x-axis label formatting depending on whether there are labels
  if (nrow(Highbarcodes) > 0) {
    bubble.plot <- bubble.plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0.5,
          colour = HighbarcodesLabel$Color,
          size = 5
        )
      )
  }

  bubble.plot <- bubble.plot +
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
