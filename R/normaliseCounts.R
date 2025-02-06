#' @title
#' Normalise barcode counts
#'
#' @description
#' Normalise barcode counts within a DGEList object
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param method method of normalization. One of `CPM`, `TMM`, `TMMwsp`, `upperquartile` or `RLE` (string). See `edgeR::calcNormFactors()`.
#' @param threshold Threshold to apply to counts before normalisation (integer). Default = `0`.
#'
#' @return Returns a DGEList object containing normalised barcode counts
#' @export
#'
#' @examples
#' data(test.dge)
#' test.dge.normalized <- normaliseCounts(test.dge, method = "CPM", threshold = 0)
#' head(test.dge.normalized)
#'

normaliseCounts <-
  function(dgeObject,
           method = "CPM",
           threshold = 0) {
    inputChecks(dgeObject)

    # check method
    if (!method %in% c("CPM", "RLE", "TMM", "upperquartile", "TMMwsp")) {
      stop("method argument must be one of CPM, RLE, TMM, upperquartile or TMMwsp")
    }

    # check threshold
    if (threshold < 0) {
      stop("threshold value should be postitive")
    }

    # threshold counts before normalisation
    if (threshold > 0) {
      message(paste("Applying raw count threshold =", threshold))
    }

    keeprows = rowSums(dgeObject$counts) >= as.numeric(threshold)
    dgeObject <- dgeObject[keeprows, ]

    # default is raw CPM
    if (method == "CPM") {
      norm.counts <- edgeR::cpm(dgeObject)
    } else {
      norm.counts <-
        cpm(edgeR::calcNormFactors(dgeObject$counts, method = method))
    }

    dgeObject.norm <- dgeObject
    dgeObject.norm$counts <- norm.counts

    return(dgeObject.norm)
  }
