#' getDominantBarcodes
#'
#' get names of barcodes over a given percentage frequency per sample
#'
#' @title
#' Get dominant barcodes per sample
#'
#' @description
#' Takes a DGEList or dataframe of barcode counts,
#' computes the percentage abundance of each barcode within each sample
#' then returns a list of barcodes meeting a threshold abundance per sample
#'
#' @param counts.obj DGEList object or dataframe containing raw / normalised counts of barcodes
#' @param pct.thresh percentage abundance threshold to call dominant barcodes
#'
#' @return Returns a named list containing vectors of dominant barcodes per sample
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' getDominantBarcodes(test.dge, pct.thresh = 5)

getDominantBarcodes <- function(counts.obj, pct.thresh = 5){

  # get counts from DGEList object if given
  if(class(counts.obj) == "DGEList"){
    counts.obj <- edgeR::cpm(counts.obj$counts)
  }

  # transform CPM into percentage within sample
  barcodes.proportional <- as.data.frame(counts.obj)
  barcodes.proportional <- sweep(barcodes.proportional,2,colSums(barcodes.proportional),`/`) * 100

  samps <- colnames(barcodes.proportional)

  # select all barcodes above threshold per sample
  dominant.barcodes <- lapply(samps, function(x){
    sample.df <- barcodes.proportional[,as.character(x), drop = F]
    sample.df <- sample.df[order(sample.df[,1], decreasing = T),,drop = F]
    dominant <- which(sample.df >= pct.thresh)
    dominant.df <- sample.df[dominant,,drop = F]
    samp.dominant <- rownames(dominant.df)
  })
  names(dominant.barcodes) <- samps

  return(dominant.barcodes)
}
