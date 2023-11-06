#' plotBarcodeHistogram
#'
#' Generate stacked barcode plots showing proportion from raw count object.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param sample Optional, name of sample to order all barcodes against.
#' @param topN top n barcodes in sample to color, all other barcodes are shown in grey (integer). Default = `10`.
#' @param title Plot title (string). Default = `Barcode Histogram`.
#'
#' @return Returns a barcode histogram plot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' data(test.dge)
#' plotBarcodeHistogram(test.dge, sample = "T0-1")

plotBarcodeHistogram <-
  function(dgeObject,
           sample = NULL,
           topN = 10,
           name = "Barcode Histogram") {

    if (methods::is(dgeObject)[1] != "DGEList") {
      stop("Please supply a valid DGEList object as input")
    }
    counts <- dgeObject$counts

    # convert counts to a dataframe of proportions
    barcodes.proportional <-
      as.data.frame(sweep(counts, 2, colSums(counts), `/`) * 100)

    if (!is.null(sample)) {
      if (length(which(
        as.character(sample) %in% colnames(barcodes.proportional)
      )) != 1) {
        stop("sample not found")
      } else {
        barcodes.proportional <-
          barcodes.proportional[order(barcodes.proportional[, sample], decreasing = T), , drop = F]
      }
    }

    top.bc <- utils::head(rowSums(barcodes.proportional), n = topN)

    # set colors
    top.colors <- grDevices::rainbow(length(names(top.bc)))
    all.colors <-
      c(top.colors, rep("grey80", length(rownames(
        barcodes.proportional
      )) - length(top.colors)))
    barcodes.proportional$color <- all.colors

    #par(mar=c(3,10,1.5,2) +.1)
    graphics::barplot(
      as.matrix(barcodes.proportional[, 1:length(names(barcodes.proportional)) -
                                        1]),
      beside = F,
      horiz = T,
      border = T,
      col = barcodes.proportional$color,
      names.arg = colnames(barcodes.proportional)[1:length(colnames(barcodes.proportional)) -
                                                    1],
      las = 1,
      cex.names = .5,
      cex.axis = 0.6,
      xlab = "Barcode proportion (%)",
      main = name
    )

  }
