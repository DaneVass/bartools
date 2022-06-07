#' plotBarcodePCA
#'
#' Plot the first two principal components of the barcode dataset
#'
#' @title
#' Barcode PCA plot
#'
#' @description
#' Takes an edgeR DGEList object contining barcode data or a dataframe of the barcode cpm counts and performs PCA.
#' Then, plots the first two dimensions.
#' Modified from DESeq2 package plotPCA function
#'
#' @param object DGEList or dataframe containing normalised barcode counts
#' @param intgroup vector of grouping variables of interest
#' @param col color group for plot, must be one of intgroup
#' @param ntop number of top most variable genes to be used in PCA calculation
#' @param returnData Logical. return a data.frame of PCA calculation?
#' @param batch metadata category indicating source of batch effects. Will be used in limma removeBatchEffect function prior to PCA
#'
#' @return Returns a plot of the first two principal components in the dataset
#' @export
#'
#' @examples
#' data(test.dge)
#' plotBarcodePCA(test.dge, intgroup = "group", ntop = 500, returnData = FALSE, batch = NULL)

plotBarcodePCA <- function(object, intgroup = "condition", col = "group", ntop = 500, returnData = FALSE, batch = NULL){

  # modified from the DEseq2 plotPCA function - LGPL license
  # https://github.com/Bioconductor-mirror/DESeq2/blob/0d7983c345bfc576725ef89addcb49c6e14ef83d/R/plots.R

  # get data from DGEList or dataframe
  if (class(object) == "DGEList"){
    dge <- edgeR::calcNormFactors(object, method = "TMM")
    data <- dge$counts * edgeR::getOffset(dge)
  } else {
    data <- as.data.frame(edgeR::cpm(object))
  }

  # get batch corrected dataset using limma removeBatchEffect
  if (!is.null(batch)){
    data <- limma::removeBatchEffect(dge, batch = batch)
  }

  # calculate the variance for each gene based on corrected dataset
  rv <- matrixStats::rowVars(data)

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # Correct for batch effects and perform a PCA on the data in assay(x) for the selected genes
  pca <- stats::prcomp(t(data[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(object$samples))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(object$samples[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    object$samples[[intgroup]]
  }

  # assembly the data for the plot
  PCAdata <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  attr(PCAdata, "percentVar") <- percentVar[1:2]

  # return the data frame if required otherwise return the plot
  if (returnData) {
    return(PCAdata)
  } else {
    percentVar <- round(100 * attr(PCAdata, "percentVar"))
    p <- ggplot2::ggplot(PCAdata, ggplot2::aes(x = PCAdata$PC1, y = PCAdata$PC2, color=PCAdata[,col], alpha=I(0.8))) +
      ggplot2::geom_point(size=3) +
      ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggplot2::theme_bw()

  }
  return(p)
}




