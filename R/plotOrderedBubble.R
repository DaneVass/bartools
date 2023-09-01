#' plotOrderedBubble
#'
#' Generate ordered proportional bubbleplots from raw count object with barcodes labelled above a specified threshold
#'
#' @param counts dataframe containing raw counts of barcodes. Expects barcodes as row names and samples as columns. Alternatively, DGE object.
#' @param proportionCutoff barcodes represented at a percentage within any sample above this threshold will be labelled
#' @param title desired plot title
#' @param orderSample name of sample to order by
#' @param samples Dataframe containing sample metadata sheet with samples as row names.
#' @param group metadata field to annotate samples on bubble plot
#' @param displaySamples vector of samples to display - keep the order of vector
#' @param displayBarcodes vector of barcodes to display
#' @param colorDominant only color clones with frequency above proportionCutoff. Others colored grey
#' @param filterCutoff barcodes below this threshold in orderSample will be filtered in all samples
#' @param legend Boolean, whether to print a legend of bubble sizes
#'
#' @return Returns a bubbleplot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' plotOrderedBubble(test.dge$counts, orderSample = "T0-1")

plotOrderedBubble <- function(counts,
                              title = "Proportional Bubble Plot",
                              orderSample = NULL,
                              samples = NULL,
                              group = NULL,
                              displaySamples = NULL,
                              displayBarcodes = NULL,
                              proportionCutoff = 10,
                              colorDominant = FALSE,
                              filterCutoff = NULL,
                              legend = TRUE) {
  ###### check inputs ##########
  if (methods::is(counts)[1] == "DGEList") {
    # if no samplesheet is provided, extract it from DGE object
    if (is.null(samples)) {
      samples <- as.data.frame(counts$samples)
    }
    counts <- as.data.frame(counts$counts)
  }
  else {
    counts <- as.data.frame(counts)
  }
  # check parameters

  if (is.null(orderSample)) {
    stop("Please provide sample to order by")
  }

  if (!is.null(group) & is.null(samples)) {
    stop("If grouping samples, samples must be provided or present in DGE object")
  }

  if (!is.null(group) & !all(group %in% colnames(samples))) {
    stop("Group must be column in samples")
  }

  if (!orderSample %in% colnames(counts)) {
    stop("Sample to order by is not present in counts data")
  }

  if (!is.null(displaySamples)) {
    # check that sample names are in count object
    missingSamples <-
      setdiff(displaySamples, colnames(counts))
    if (length(missingSamples) > 0) {
      stop(paste(
        "Samples",
        paste(missingSamples, collapse = ", "),
        "not found in count object"
      ))
    }
  }

  if (!is.null(displayBarcodes)) {
    # check that sample names are in count object
    missingBarcodes <-
      setdiff(displayBarcodes, rownames(counts))
    if (length(missingBarcodes) > 0) {
      stop(paste(
        "Barcodes",
        cat(missingBarcodes),
        "not found in count object"
      ))
    }
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
  # high prop barcodes
  # TODO
  Highbarcodes <-
    dplyr::filter_all(barcodes.proportional[, 1:(ncol(barcodes.proportional) - 3)],
                      dplyr::any_vars(. > proportionCutoff))

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

  # Assign diverse colours to high frequency barcodes
  SelColors <- scales::hue_pal()(nrow(Highbarcodes))
  i = 1
  for (bc in rownames(Highbarcodes)) {
    barcodes.proportional[bc,]$Color <- SelColors[i]
    i <- i + 1
  }
  HighbarcodesLabel <-
    barcodes.proportional[rownames(Highbarcodes),]
  HighbarcodesLabel <-
    HighbarcodesLabel[HighbarcodesLabel$Position > 0,]


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
                  x = "",
                  title = title) +
    ggplot2::scale_size_continuous(
      name = "Barcode Proportion (%)",
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
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        colour = HighbarcodesLabel$Color,
        size = 6
      ),
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
