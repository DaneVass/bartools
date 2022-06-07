#' Normalise barcode counts
#'
#' Normalise barcode counts within a DGEList object
#'
#' @param dge DGEList object containing raw barcode counts to be normalised
#' @param method method of normalisation
#' @param threshold Threshold to apply to counts before normalisation
#'
#' @return Returns a data.frame object containing normalised barcode counts
#' @export
#'
#' @examples
#' data(test.dge)
#' normaliseCounts(test.dge, method = "cpm", threshold = 0)
#'

normaliseCounts <- function(dge, method = "cpm", threshold = 0){

  # threshold counts before normalisation
  keeprows = rowSums(dge$counts) >= as.numeric(threshold)
  dge <- dge[keeprows,]

  if(method == "cpm"){
    norm.counts <- as.data.frame(edgeR::cpm(dge))
  }

  return(norm.counts)
}
