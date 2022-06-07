#' calcPercentileBarcodes
#'
#' Calculate barcodes comprising the Nth percentile for each sample & generate cumulative sum plots
#'
#' @param counts.obj DGEList object containing raw or normalised barcode counts
#' @param percentile desired percentile value
#'
#' @return Returns a list object containing nth percentile tables of barcode counts per sample
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' calcPercentileBarcodes(test.dge, percentile = .95)

calcPercentileBarcodes <- function(counts.obj, percentile = .95){

  counts.obj <- as.data.frame(counts.obj$counts)
  samples <- colnames(counts.obj)

  # setup empty dataframes to collect cumulative sum data
  TotalNumPercentile <- data.frame()
  BarcodesPercentile <- vector(mode = "list", length = length(samples))
  BarcodesCountPercentile <- vector(mode = "list", length = length(samples))

  # for each sample calculate the nth percentile barcodes.
  for (i in 1:ncol(counts.obj)){

    # sort dataset
    ordered <- order(counts.obj[,i], decreasing = T)
    sorted <- counts.obj[ordered,i,drop = F]
    sorted <- sorted[sorted > 0,, drop = F]

    colsum <- sum(sorted)
    percentile.cutoff <- sum(sorted)*percentile

    # generate cumulative sum of sorted dataset
    cumsum <- cumsum(sorted)

    # find number of barcodes that make up Nth percentile
    len <- length(which(cumsum <= percentile.cutoff))
    print(paste(samples[i], " ", percentile, "th percentile: ", len))

    # fill Nth percentile data frame for plotting below
    d <- data.frame(Sample=factor(samples[i]),NumBarcodes=len)
    TotalNumPercentile <- rbind(TotalNumPercentile, d)

    # extract barcodes in nth percentile
    topBC <- rownames(sorted)[which(cumsum <= percentile.cutoff)]
    BarcodesPercentile[[i]] <- topBC
    names(BarcodesPercentile) <- samples

    # table of barcodes with numbers in nth percentile
    BarcodesCountPercentile[[i]] <- sorted[1:length(topBC),,drop=F]
    names(BarcodesCountPercentile) <- samples
  }
  return (list("NumBarcodes" = TotalNumPercentile,
               "TopBarcodes" = BarcodesPercentile,
               "TopBarcodeCounts" = BarcodesCountPercentile))

}


