#' plotBarcodePCA
#'
#' Plot the first two principal components of the barcode dataset
#'
#' @title
#' Barcode PCA plot
#'
#' @description
#' Takes a DGEList object contining barcode cpm counts and performs PCA.
#' Then, plots the first two dimensions.
#' Modified from DESeq2 package plotPCA function
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param groups Optional, one or multiple column names in sample metadata to color samples by (string or vector of strings).
#' @param ntop number of top most variable barcodes to be used in PCA calculation (integer). Default = `500`.
#' @param returnData Return a data.frame of PCA calculation instead of plot (boolean). Default = `FALSE`.
#' @param dimensions Principle components to plot (vector of integers of length 2). Default = `c(1, 2)`.
#' @param batch Optional, metadata category indicating source of batch effects. Will be used in limma removeBatchEffect function prior to PCA.
#'
#' @return Returns a plot of the first two principal components in the dataset
#' @export
#'
#' @examples
#' plotBarcodePCA(test.dge, groups = "Treatment", ntop = 500, returnData = FALSE, batch = NULL)

plotBarcodePCA <-
  function(dgeObject,
           groups = NULL,
           ntop = 500,
           returnData = FALSE,
           pcs = c(1, 2),
           batch = NULL) {
    inputChecks(dgeObject, groups = groups)
    # modified from the DEseq2 plotPCA function - LGPL license
    # https://github.com/Bioconductor-mirror/DESeq2/blob/0d7983c345bfc576725ef89addcb49c6e14ef83d/R/plots.R

    # get data from DGEList
    dge <- edgeR::calcNormFactors(dgeObject, method = "TMM")

    # get batch corrected dataset using limma removeBatchEffect
    if (!is.null(batch)) {
      data <- limma::removeBatchEffect(dge, batch = batch)
    } else {
      data <- dge$counts * edgeR::getOffset(dge)
    }

    # calculate the variance for each gene based on corrected dataset
    rv <- apply(data, 1, stats::var)

    # select the ntop barcodes by variance
    selected_barcodes <-
      order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

    # Correct for batch effects and perform a PCA on the data in assay(x) for the selected barcodes
    pca <- stats::prcomp(t(data[selected_barcodes,]))

    # the contribution to the total variance for each component
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)

    intgroup.df <-
      as.data.frame(dgeObject$samples[, groups, drop = FALSE])

    # add the intgroup factors together to create a new grouping factor
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))

    # assembly the data for the plot
    PCAdata <-data.frame(pca$x, group = group,
               intgroup.df,
               name = colnames(dgeObject))

    attr(PCAdata, "percentVar") <- percentVar

    # return the data frame if required otherwise create the plot
    if (returnData) {
      return(PCAdata)
    }

    percentVar <- round(100 * attr(PCAdata, "percentVar"))
    p <-
      ggplot2::ggplot(PCAdata,
                      ggplot2::aes_string(
                        x = paste0("PC", pcs[1]),
                        y = paste0("PC", pcs[2]),
                        color = "group"
                      )) +
      ggplot2::geom_point(size = 3, alpha=0.7) +
      ggplot2::xlab(paste0("PC", pcs[1], ": ", percentVar[pcs[1]], "% variance")) +
      ggplot2::ylab(paste0("PC", pcs[2], ": ", percentVar[pcs[2]], "% variance")) +
      ggplot2::theme_bw() +
      ggplot2::guides(color = guide_legend(title = paste(groups, collapse = " : ")))

    return(p)
  }
