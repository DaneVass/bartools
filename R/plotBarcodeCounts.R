#' Plot total counts per barcode in library
#'
#' Simple plot of total read counts per barcode in library
#'
#' @param counts data.frame of barcode count x sample
#' @param order Logical. Order the dataset be decreasing abundance
#' @param log10 Logical. log10 transform the data. Adds a pseudocount of 1
#'
#' @return Returns a plot of the read counts per barcode (row) in a data frame
#' @import ggplot2
#' @export
#'
#' @examples
#' data(test.counts)
#' plotBarcodeCounts(test.counts)
#' plotBarcodeCounts(test.counts, order = TRUE)
#' plotBarcodeCounts(test.counts, order = TRUE, log10 = TRUE)

plotBarcodeCounts <- function(counts, order = F, log10 = F){
  rowsums <- rowSums(counts)

  if(log10){
    rowsums <- log10(rowsums+1)
  }

  if(order){
    ordered <- sort(rowsums, decreasing = T)
    graphics::barplot(ordered,
            las=2,
            main="Total counts per barcode",
            axisnames = F,
            cex.axis=0.8,
            xlab = "Barcode - Descending order by total read count across samples", ylab = "Barcode total read count")
  } else {
    graphics::barplot(rowsums,
            las=2,
            main="Total counts per barcode",
            axisnames = F,
            cex.axis=0.8,
            xlab = "Barcode - Descending order by frequency in reference library", ylab = "Barcode total read count")

  }
}
