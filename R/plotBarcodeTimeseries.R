#' plotBarcodeTimeseries
#'
#' Generate proportional timeseries plot from raw / normalised barcode count object.
#'
#' @param counts.obj dataframe containing raw counts of barcodes
#' @param name desired plot title
#' @param seed RNG seed
#' @param top number of top barcodes per sample to plot

#' @return Returns a bubbleplot of barcodes represented by proportion of total pool

#' @export

#' @examples
#' data(test.counts)
#' plotBarcodeTimeseries(test.counts[,1:4], name = "Proportional Timeseries Plot", seed = 5)

plotBarcodeTimeseries <- function(counts.obj, name = "", seed = 5, top = 50){
    barcodes.proportional <- as.data.frame(counts.obj)
    barcodes.proportional <- sweep(barcodes.proportional,2,colSums(barcodes.proportional),`/`) * 100
    barcodes.proportional$barcode <- rownames(barcodes.proportional)

    # take top n barcodes per sample
    top.bc <- lapply(colnames(barcodes.proportional)[-length(colnames(barcodes.proportional))],
                     function(x){
                       dat <- dplyr::select(barcodes.proportional, `barcode`, as.character(x))
                       dat <- dat[order(dat[,2], decreasing = T),][1:top,]
                       return(dat$barcode)
    })

    names(top.bc) <- colnames(barcodes.proportional)[-length(colnames(barcodes.proportional))]

    # get unique list of barcodes
    top.bc <- unique(unlist(top.bc))

    # subset counts data to keep only top bc
    barcodes.proportional <- barcodes.proportional[top.bc,]
    barcodes.proportional.melted <- reshape2::melt(barcodes.proportional)

    colnames(barcodes.proportional.melted) <- c("Barcode", "Sample", "Proportion")
    barcodes.proportional.melted$Barcode <- as.factor(barcodes.proportional.melted$Barcode)
    barcodes.proportional.melted$Proportion <- as.numeric(barcodes.proportional.melted$Proportion)

    timepoints <- unique(barcodes.proportional.melted$Sample)

    colors <- scales::hue_pal()(30)
    set.seed(seed) # set custom seed to get same color order every time
    colors <- sample(colors, length(rownames(barcodes.proportional.melted)), replace = TRUE)
    barcodes.proportional.melted$color <- colors

    # inspired by genBaRcode package
    timeseries.plot <- ggplot2::ggplot(barcodes.proportional.melted,
                              ggplot2::aes_string(group = "Barcode",fill = "Barcode", x = "Sample", y = "Proportion")) +
        ggplot2::theme_bw() +
        ggplot2::geom_area(alpha=0.9) +
        ggplot2::scale_fill_manual(values = colors) +
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_x_discrete(breaks = timepoints, labels = timepoints, limits = timepoints) +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       panel.grid.major.y = element_blank(),
                       panel.grid.minor.y = element_blank()) +
        ggplot2::ggtitle(label = paste("Proportional Timeseries Plot:", name))

    return(timeseries.plot)
}
