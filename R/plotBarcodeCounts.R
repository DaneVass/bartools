#' Plot total counts per barcode in library
#'
#' Simple plot of total read counts per barcode across samples.
#'
#' @param dgeObject DGE object with barcode counts
#' @param order Order the dataset by decreasing abundance (boolean). Default = `FALSE`.
#' @param log10 log10 transform the data. Adds a pseudocount of 1 (boolean). Default = `FALSE`.
#'
#' @return Returns a plot of the read counts per barcode (row) in a data frame
#' @import ggplot2
#' @export
#'
#' @examples
#' data(test.dge)
#' plotBarcodeCounts(test.dge)
#' plotBarcodeCounts(test.dge, order = TRUE, log10 = TRUE)

plotBarcodeCounts <-
  function(dgeObject,
           order = FALSE,
           log10 = FALSE) {
    inputChecks(dgeObject)
    counts <- dgeObject$counts

    # sort barcodes alpha-numerically
    counts <- counts[stringr::str_sort(rownames(counts), numeric = T), ]

    rowsums <- rowSums(counts)

    if (log10) {
      rowsums <- log10(rowsums + 1)
    }

    if (order) {
      ordered <- sort(rowsums, decreasing = T)
      graphics::barplot(
        ordered,
        las = 2,
        main = "Total counts per barcode",
        axisnames = F,
        cex.axis = 0.8,
        xlab = "Barcode - Descending order by total read count across samples",
        ylab = "Barcode total read count"
      )
    } else {
      graphics::barplot(
        rowsums,
        las = 2,
        main = "Total counts per barcode",
        axisnames = F,
        cex.axis = 0.8,
        xlab = "Barcode - Descending order by frequency in reference library",
        ylab = "Barcode total read count"
      )

    }
  }
