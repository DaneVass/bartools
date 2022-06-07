#' plotBCBoxplot
#'
#' Boxplot of selected barcodes
#'
#' @title
#' Boxplot of selected barcodes
#'
#' @description
#' Takes a dataframe of barcode counts,
#' select specific barcodes in the dataframe
#' and boxplot their counts
#'
#' @param dge.obj DGE object containing counts and sample information
#' @param barcodes vector of barcodes to be plotted
#' @param condition specific metadata condition to be plotted
#'
#' @return Returns a boxplot
#' @export
#' @examples
#' data(test.dge)
#' barcodes <- sample(rownames(test.dge), 10)
#' plotBarcodeBoxplot(test.dge, barcodes)
#' plotBarcodeBoxplot(test.dge, barcodes, condition = c("T0", "T0b"))

plotBarcodeBoxplot <- function(dge.obj, barcodes, condition = NULL){
  counts = dge.obj$counts
  sample = dge.obj$samples
  # no condition specified
  if(is.null(condition)){
    # only one barcode given
    if (length(barcodes) == 1){
      graphics::boxplot(counts[rownames(counts) %in% barcodes,], ylab = "Barcode counts", xlab = barcodes[1])
    } else{
      tab <- t(counts[rownames(counts) %in% barcodes,])
      stacked_tab <- reshape2::melt(tab)
      graphics::par(mar=c(10,5,2,2))
      graphics::boxplot(stacked_tab[,3]~ stacked_tab[,2],col = grDevices::rainbow(ncol(tab)),xlab="",ylab = "Barcode counts",las=3)
    }
  } else{
    # number of unique condition to be plotted
    L = length(unique(condition))
    graphics::par(mfrow=c(ceiling(L/3),L))
    graphics::par(mar=c(10,3,1,1))
    if (length(barcodes) > 1){
      for (i in 1:L){

        c1 = unique(condition)[i]
        s1 = rownames(sample[condition == c1,])
        tab = t(counts[rownames(counts) %in% barcodes,colnames(counts) %in% s1])
        stacked_tab <- reshape2::melt(tab)
        graphics::boxplot(stacked_tab[,3]~ stacked_tab[,2],col = grDevices::rainbow(ncol(tab)),xlab="",ylab = "Barcode counts",las=3,main=c1)
      }
    } else{
      # only one barcode
      for (i in 1:L){
        c1 = unique(condition)[i]
        s1 = rownames(sample[condition == c1,])
        tab = counts[rownames(counts) %in% barcodes,colnames(counts) %in% s1]
        graphics::boxplot(tab,xlab="",ylab = "Barcode counts",las=3,main=c1)
      }
    }
  }
}
