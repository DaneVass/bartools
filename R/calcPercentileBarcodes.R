#' calcPercentileBarcodes
#'
#' Calculate barcodes comprising the top Nth percentile for each sample.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param percentile Desired percentile value (decimal). Default = `0.95` (i.e. 95th percentile).
#'
#' @return Returns a list object containing Nth percentile tables of barcode counts per sample
#' - TotalNumPercentile: Table with the number of barcodes in the top Nth percentile per sample
#' - TopBarcodes: List of barcodes in the top Nth percentile per sample
#' - TopBarcodeCounts: List of tables with barcodes and counts in the top Nth percentile per sample
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' calcPercentileBarcodes(test.dge, percentile = 0.95)

calcPercentileBarcodes <- function(dgeObject, percentile = 0.95) {
  ###### check inputs ##########
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
    # order barcodes by decreasing counts
    ordered <- order(counts[, i], decreasing = T)
    sorted <- counts[ordered, i, drop = F]
    # remove barcodes with 0 counts
    sorted <- sorted[sorted > 0, , drop = F]

    percentile.cutoff <- sum(sorted) * percentile

    # calculate cummulative sum of barcode counts
    cumsum <- cumsum(sorted)

    # find number of barcodes that make up Nth percentile
    # check how many barcodes have a cummulative sum smaller or equal to cutoff
    len <- length(which(cumsum <= percentile.cutoff))
    # the number of barcodes should be inclusive,
    # i.e. at least 1 barcode makes up the top x percentile
    # unless the barcode is exactly at the cutoff, add 1
    if (len == 0) {
      len <- len + 1
    } else if (cumsum[len, ] != percentile.cutoff &
               len <= nrow(sorted)) {
      len <- len + 1
    }

    # fill Nth percentile data frame for plotting below
    d <- data.frame(Sample = factor(samples[i]), NumBarcodes = len)
    TotalNumPercentile <- rbind(TotalNumPercentile, d)

    # extract barcodes in nth percentile
    topBC <- rownames(sorted)[seq(len)]
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
