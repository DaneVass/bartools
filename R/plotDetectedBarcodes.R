#' plotDetectedBarcodes
#'
#' Plot the total number of barcodes detected in a sample
#'
#' @param counts DGEList or dataframe containing raw or normalised barcode counts
#' @param percentile desired percentile value. 95th percentile by default
#' @param plot Logical. plot data instead of returning counts table.
#'
#' @return Returns a histogram plot of the number of detected barcodes per sample.
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' plotDetectedBarcodes(test.dge, percentile = .95)
#' plotDetectedBarcodes(test.dge, plot = FALSE)

plotDetectedBarcodes <- function(counts, percentile = .95, plot = T){

  if (class(counts) == "DGEList"){
    counts.obj <- as.data.frame(counts$counts)
  } else {
    counts.obj <- as.data.frame(counts)
  }

  dim(counts.obj)
  samples <- colnames(counts.obj)

  # setup empty dataframe and vector to collect cumulative sum data
  percentile.df <- data.frame()
  barcodes <- c()

  for (i in samples){

    # error handling for colors field
    if(i == 'color'){
      message("skipping color column")
    }

    else {

      # sort dataset
      order <- order(counts.obj[,i], decreasing = T)
      sorted <- counts.obj[order,as.character(i), drop = F]
      sorted <- sorted[sorted > 0, , drop = F]

      # calculate percentile cutoff
      colsum <- sum(sorted)
      percentile.cutoff <- sum(sorted) * percentile

      # generate cumulative sum of sorted dataset
      cumsum <- cumsum(sorted)
      stopifnot(max(cumsum) == colsum)

      # find number of barcodes that make up Nth percentile
      len <- length(which(cumsum <= percentile.cutoff))

      # fill Nth percentile data frame for plotting below
      d <- data.frame(Sample=factor(i),Barcodes=len)
      percentile.df <- rbind(percentile.df, d)

      # extract barcodes in Nth percentile
      rows <- rownames(sorted)[which(cumsum <= percentile.cutoff)]
      length(rows)
      sorted.top <- counts.obj[rows,i,drop = F]
      top.bc.sum <- sum(sorted.top[,1])
      sorted.top$percentage <- (sorted.top[,1]/top.bc.sum)*100
    }
  }

  if (plot){
    # return a plot of the data
    p <- ggplot2::ggplot(data=percentile.df, ggplot2::aes(x = percentile.df$Sample,
                                                          y = percentile.df$Barcodes)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::guides(fill="none") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=-.5, hjust = 1)) +
      ggplot2::ggtitle("Number of detected barcodes per sample") +
      ggplot2::ylim(0,max(percentile.df$Barcodes+25)) +
      ggplot2::geom_text(ggplot2::aes(label = percentile.df$Barcodes, y = percentile.df$Barcodes + 0.05), position = ggplot2::position_dodge(0.9), vjust = -0.2)

    return(p)

  } else {
    return(percentile.df)
  }
}
