#' proportionalBubbleplot
#'
#' Generate proportional bubbleplots from raw count object with barcodes labelled above a specified threshold
#'
#' @param counts.obj dataframe containing raw counts of barcodes. assumes barcodes as rownames.
#' @param labels Boolean. print barcode labels?
#' @param proportion.cutoff barcodes represented at a percentage within any sample above this threshold will be labelled
#' @param name desired plot title
#'
#' @return Returns a bubbleplot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' data(test.counts)
#' proportionalBubbleplot(test.counts, name = "Proportional Bubble Plot", proportion.cutoff = 10)

proportionalBubbleplot <- function(counts.obj, labels = TRUE, proportion.cutoff = 10, name = "Proportional Bubble Plot"){
  # transform CPM into percentage within sample
  barcodes.proportional <- as.data.frame(counts.obj)
  barcodes.proportional <- sweep(barcodes.proportional,2,colSums(barcodes.proportional),`/`) * 100

  # give each barcode a specific color
  colors <- scales::hue_pal()(30)
  colors <- sample(colors, length(rownames(barcodes.proportional)), replace = TRUE)
  names(colors) <- rownames(barcodes.proportional)
  barcodes.proportional$Color <- colors

  # maintain Rank information of barcodes. Plot in ascending rank order from original barcode library
  barcodes.proportional$Position <- as.factor(seq(1,length(rownames(counts.obj))))
  barcodes.proportional$Barcode <- rownames(barcodes.proportional)

  # melt data frame and rename columns correctly
  barcodes.proportional.melted <- reshape2::melt(barcodes.proportional, id.vars = c("Color", "Position", "Barcode"))
  colnames(barcodes.proportional.melted) <- c("Color", "Position", "Barcode", "Sample", "Proportion")
  # convert variables to correct form
  barcodes.proportional.melted$Barcode <- as.factor(barcodes.proportional.melted$Barcode)
  barcodes.proportional.melted$Sample <- as.factor(barcodes.proportional.melted$Sample)
  barcodes.proportional.melted$Proportion <- as.numeric(barcodes.proportional.melted$Proportion)

  if(isTRUE(labels)){
    # identify high proportion barcodes for labelling
    Highbarcodes <- barcodes.proportional.melted[barcodes.proportional.melted$Proportion > proportion.cutoff,]
    # only take unique barcode labels
    HighbarcodeOrdered <- unique(Highbarcodes[order(Highbarcodes$Barcode),1:3])

    # generate bubbleplot
    bubble.plot <- ggplot2::ggplot(barcodes.proportional.melted,
                                   ggplot2::aes(x=barcodes.proportional.melted$Position,
                                                y=barcodes.proportional.melted$Sample,
                                                size = barcodes.proportional.melted$Proportion,
                                                color = barcodes.proportional.melted$Color)) +
      ggplot2::geom_point(stat = "identity", alpha = 0.6, shape = 16) +
      ggplot2::scale_color_identity() +
      ggplot2::labs(y = "Condition", x = "Barcode", title = name) +
      ggplot2::scale_size_continuous(range = c(0.1, 10)) +
      ggplot2::scale_x_discrete(breaks = HighbarcodeOrdered$Position, labels = HighbarcodeOrdered$Barcode) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
            axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, colour = HighbarcodeOrdered$Color, size = 5),
            axis.text.y = ggplot2::element_text(size = 6),
            plot.title = ggplot2::element_text(size=8),
            axis.title = ggplot2::element_text(size= 6))
  } else {
    # generate bubbleplot
    bubble.plot <- ggplot2::ggplot(barcodes.proportional.melted,
                                   ggplot2::aes(x=barcodes.proportional.melted$Position,
                                                y=barcodes.proportional.melted$Sample,
                                                size = barcodes.proportional.melted$Proportion,
                                                color = barcodes.proportional.melted$Color)) +
      ggplot2::geom_point(stat = "identity", alpha = 0.6, shape = 16) +
      ggplot2::scale_color_identity() +
      ggplot2::labs(y = "Condition", x = "Barcode", title = name) +
      ggplot2::scale_size_continuous(range = c(0.1, 10)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = 6),
            plot.title = ggplot2::element_text(size=8),
            axis.title = ggplot2::element_text(size= 6))
  }
  print(bubble.plot)
}
