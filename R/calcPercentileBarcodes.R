#' calcPercentileBarcodes
#'
#' Calculate barcodes comprising the Nth percentile for each sample & generate cumulative sum plots
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param percentile Desired percentile value (decimal). Default = `0.95` (i.e. 95 percent).
#'
#' @return Returns a list object containing Nth percentile tables of barcode counts per sample
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' calcPercentileBarcodes(test.dge, percentile = 0.95)

calcPercentileBarcodes <- function(dgeObject, percentile = 0.95) {
  inputChecks(dgeObject)

  counts <- as.data.frame(dgeObject$counts)

  samples <- colnames(counts)

  # setup empty dataframes to collect cumulative sum data
  TotalNumPercentile <- data.frame()
  BarcodesPercentile <-
    vector(mode = "list", length = length(samples))
  BarcodesCountPercentile <-
    vector(mode = "list", length = length(samples))

  # for each sample calculate the nth percentile barcodes.
  for (i in 1:ncol(counts)) {
    # sort dataset
    ordered <- order(counts[, i], decreasing = T)
    sorted <- counts[ordered, i, drop = F]
    sorted <- sorted[sorted > 0, , drop = F]

    percentile.cutoff <- sum(sorted) * percentile

    # generate cumulative sum of sorted dataset
    cumsum <- cumsum(sorted)

    # find number of barcodes that make up Nth percentile
    len <- length(which(cumsum <= percentile.cutoff))

    # fill Nth percentile data frame for plotting below
    d <- data.frame(Sample = factor(samples[i]), NumBarcodes = len)
    TotalNumPercentile <- rbind(TotalNumPercentile, d)

    # extract barcodes in nth percentile
    topBC <- rownames(sorted)[which(cumsum <= percentile.cutoff)]
    BarcodesPercentile[[i]] <- topBC
    names(BarcodesPercentile) <- samples

    # table of barcodes with numbers in nth percentile
    BarcodesCountPercentile[[i]] <- sorted[1:length(topBC), , drop = F]
    names(BarcodesCountPercentile) <- samples
  }
  return(
    list(
      "NumBarcodes" = TotalNumPercentile,
      "TopBarcodes" = BarcodesPercentile,
      "TopBarcodeCounts" = BarcodesCountPercentile
    )
  )

}
