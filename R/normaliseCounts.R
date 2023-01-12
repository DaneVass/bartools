#' Normalise barcode counts
#'
#' Normalise barcode counts within a DGEList object
#'
#' @param dge DGEList object containing raw barcode counts to be normalised
#' @param method method of normalization. One of cpm, tmm, tmmwsp, upperquartile or rle. see edgeR::calcNormFactors()
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

  # check obj
  if (class(dge)[1] != "DGEList"){
    stop("Please provide a DGEList object for normalisation")
  }
  
  # check method
  if(!method %in% c("CPM", "RLE", "TMM", "upperquartile", "TMMwsp")){
    stop("method argument must be one of CPM, RLE, TMM, upperquartile or TMMwsp")
  }
  
  # check threshold
  if(threshold < 0){
    stop("threshold value shupld be postitive")
  }
  
  # threshold counts before normalisation
  if (threshold > 0){
    message(paste("Applying raw count threshold =", threshold))
  }
  
  keeprows = rowSums(dge$counts) >= as.numeric(threshold)
  dge <- dge[keeprows,]

  if(method == "CPM"){
    norm.counts <- as.data.frame(edgeR::cpm(dge))
  } else {
    norm.counts <- as.data.frame(dge$counts %*% diag(edgeR::calcNormFactors(dge$counts, method = method)))
  }

  return(norm.counts)
}
