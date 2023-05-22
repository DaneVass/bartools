#' plotOrderedBubble
#'
#' Generate ordered proportional bubbleplots from raw count object with barcodes labelled above a specified threshold
#'
#' @param counts.obj dataframe containing raw counts of barcodes. assumes barcodes as rownames.
#' @param labels logical. Print barcode labels
#' @param proportion.cutoff barcodes represented at a percentage within any sample above this threshold will be labelled
#' @param name desired plot title
#' @param orderSample name of sample to order by
#' @param colorDominant only color clones with frequency above proportion.cutoff. Others colored grey
#'
#' @return Returns a bubbleplot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' plotOrderedBubble(test.dge$counts, orderSample = "T0-1")

plotOrderedBubble <- function(counts.obj,
                              labels = T,
                              name = "Proportional Bubble Plot",
                              orderSample = NULL,
                              proportion.cutoff = 10,
                              colorDominant = FALSE) {
  if (is.null(orderSample)) {
    stop("Please provide sample to order by")
  }

  # transform CPM into percentage within sample
  barcodes.proportional <- as.data.frame(counts.obj)
  barcodes.proportional <-
    sweep(barcodes.proportional,
          2,
          colSums(barcodes.proportional),
          `/`) * 100

  # Order by selected sample. Plot in ascending rank order from original barcode library
  barcodes.proportional$Position <-
    barcodes.proportional[, orderSample]
  barcodes.proportional$Barcode <- rownames(barcodes.proportional)

  if (colorDominant) {
    # make all barcodes grey and only color those that are above threshold cutoff
    colors <- "grey90"
    colors <-
      sample(colors, length(rownames(barcodes.proportional)), replace = TRUE)
    names(colors) <- rownames(barcodes.proportional)
    barcodes.proportional$Color <- colors

    # Assign diverse colours to high frequency barcodes
    Highbarcodes <-
      dplyr::filter_all(barcodes.proportional[, 1:(ncol(barcodes.proportional) -
                                                     3)],
                        dplyr::any_vars(. > proportion.cutoff))
    SelColors <- scales::hue_pal()(nrow(Highbarcodes))
    i = 1
    for (bc in rownames(Highbarcodes)) {
      barcodes.proportional[bc, ]$Color <- SelColors[i]
      i <- i + 1
    }
  } else {
    # give each barcode a specific color
    colors <- scales::hue_pal()(30)
    colors <-
      sample(colors, length(rownames(barcodes.proportional)), replace = TRUE)
    names(colors) <- rownames(barcodes.proportional)
    barcodes.proportional$Color <- colors

    # Assign diverse colours to high frequency barcodes
    Highbarcodes <-
      dplyr::filter_all(barcodes.proportional[, 1:(ncol(barcodes.proportional) -
                                                     3)],
                        dplyr::any_vars(. > proportion.cutoff))
    SelColors <- scales::hue_pal()(nrow(Highbarcodes))
    i = 1
    for (bc in rownames(Highbarcodes)) {
      barcodes.proportional[bc, ]$Color <- SelColors[i]
      i <- i + 1
    }
  }

  HighbarcodesLabel <-
    barcodes.proportional[rownames(Highbarcodes), ]
  HighbarcodesLabel <-
    HighbarcodesLabel[HighbarcodesLabel$Position > 0, ]

  # melt data frame and rename columns correctly
  barcodes.proportional.melted <-
    reshape2::melt(barcodes.proportional,
                   id.vars = c("Color", "Position", "Barcode"))
  colnames(barcodes.proportional.melted) <-
    c("Color", "Position", "Barcode", "Sample", "Proportion")
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
    ggplot2::labs(y = "Condition", x = "Barcode Proportion", title = name) +
    ggplot2::scale_size_continuous(range = c(0.1, 10)) +
    ggplot2::scale_x_continuous(
      trans = 'log10',
      labels =
        HighbarcodesLabel$Barcode,
      breaks = HighbarcodesLabel$Position,
      sec.axis = ggplot2::sec_axis(
        ~ . * 1,
        labels = c(0.0001, 0.001, 0.01, 0.1, 1, 2, 5, 10, 20, 30, 40, 50),
        breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 2, 5, 10, 20, 30, 40, 50)
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        colour = HighbarcodesLabel$Color,
        size = 6
      ),
      axis.text.y = ggplot2::element_text(size = 6),
      plot.title = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 6)
    )

  return(bubble.plot)

  # SelColors <- data.frame(value = SelColors)
  # return(SelColors)
}