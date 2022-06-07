#' plotBarcodeHistogram
#'
#' Generate stacked barcode plots showing proportion from raw count object.
#'
#' @param counts.obj dataframe containing raw counts of barcodes
#' @param sample desired name of sample to order all barcodes against
#' @param top top n barcodes in sample to color, all other barcodes are shown in grey
#' @param name desired plot title
#'
#' @return Returns a barcode histogram plot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' data(test.counts)
#' plotBarcodeHistogram(test.counts, sample = "T0-1", top = 10, name = "Barcode Histogram")

plotBarcodeHistogram <- function(counts.obj, sample = NULL, top = 10, name = "Barcode Histogram"){

  # convert counts to a dataframe of proportions
  barcodes.proportional <- as.data.frame(sweep(counts.obj,2,colSums(counts.obj),`/`) * 100)

  if( !is.null(sample) ){
    if ( length(which(as.character(sample) %in% colnames(barcodes.proportional))) != 1 ) {
      stop("sample not found")
    } else {
      barcodes.proportional <- barcodes.proportional[order(barcodes.proportional[,sample], decreasing = T),,drop = F]
    }
  }

  top.bc <- utils::head(rowSums(barcodes.proportional), n = top)

  # set colors
  top.colors <- grDevices::rainbow(length(names(top.bc)))
  all.colors <- c(top.colors, rep("grey80", length(rownames(barcodes.proportional)) - length(top.colors)))
  barcodes.proportional$color <- all.colors

  #par(mar=c(3,10,1.5,2) +.1)
  graphics::barplot(as.matrix(barcodes.proportional[,1:length(names(barcodes.proportional))-1]),
          beside = F, horiz = T, border = T, col = barcodes.proportional$color,
          names.arg = colnames(barcodes.proportional)[1:length(colnames(barcodes.proportional))-1],
          las=1, cex.names = .5, cex.axis = 0.6, xlab = "Barcode proportion (%)", main = name)

}

