#' plotLibraryDiversity
#'
#' Lineplot of the barcode abundances in two different experimental settings
#'
#' @title
#' Lineplot of barcode abundances in two conditions
#'
#' @description
#' Takes a dataframe of barcode counts,
#' computes and plots the abundance of each barcode in the library sample
#'
#' @param barcodes dataframe containing raw counts of barcodes
#' @param samplename sample condition of interest
#' @param cutoff rowsum cutoff defining rows to keep
#' @param skew Logical. Calculate and plot the library skew ratio (90th percentile / 10th percentile)
#'
#' @return Returns a barcode library frequency plot
#' @export

# Barcode frequency distribution
plotLibraryDiversity <- function(barcodes, samplename = "Library", cutoff = 10, skew = F){

  # filter barcodes
  colnames(barcodes) <- c("Barcode", "Raw_count")
  barcodes$Raw_count <- as.numeric(barcodes$Raw_count)
  barcodes <- barcodes[which(barcodes$Raw_count >= cutoff),] # arbitrary 10 as default

  if (skew) {
    # skew ratio of top 10% to bottom 10% of guide counts
    top_10 <- stats::quantile(barcodes[,2], probs = 0.9)
    bottom_10 <- stats::quantile(barcodes[,2], probs = 0.1)
    if (top_10 != 0 && bottom_10 != 0) {
      skew_ratio = top_10/bottom_10
    } else {
      stop('Not enough perfect matches to determine skew ratio')
    }
  }

  p <- ggplot2::ggplot(barcodes, ggplot2::aes(y=barcodes$Raw_count, x=seq(1,length(rownames(barcodes))))) +
    ggplot2::geom_point(stat = "identity", show.legend = F, alpha = 0.4, size = 0.1) +
    ggplot2::scale_y_continuous(trans='log10') +
    ggplot2::theme_bw() +
    #scale_size_manual(values=c(2,2)) +
    ggplot2::geom_hline(yintercept = stats::median(barcodes$Raw_count), color = "grey40") +
    ggplot2::xlab("Barcode") +
    ggplot2::ylab("Count") +
    ggplot2::ggtitle(paste(samplename, ": Barcode frequency distribution"))

  if (skew) {
    p <- p +
      ggplot2::geom_vline(xintercept = top_10, color = "blue") +
      ggplot2::geom_vline(xintercept = bottom_10, color = "red")
  }

  return(p)
}

#' plotLibraryCumSum
#'
#' Plot the cumulative sum of the barcode abundances in a barcode library
#'
#' @title
#' Cumulative sum plot of barcode abundances in a reference library
#'
#' @description
#' Takes a dataframe of barcode counts,
#' computes the median abundance of each barcode for two specific conditions,
#' then do a line plot for both conditions
#'
#' @param barcodes dataframe containing raw counts of barcodes
#' @param samplename sample condition of interest
#' @param cutoff rowsum cutoff defining rows to keep
#'
#' @return Returns a cumulative sum plot for the desired barcode library sample
#' @export

plotLibraryCumSum <- function(barcodes, samplename = "Library", cutoff = 10){

  # filter barcodes
  colnames(barcodes) <- c("Barcode", "Raw_count")
  barcodes$Raw_count <- as.numeric(barcodes$Raw_count)
  barcodes <- barcodes[which(barcodes$Raw_count >= cutoff),] # arbitrary 10 as default

  # generate cumsum data
  sorted <- sort(barcodes$Raw_count, decreasing = T)
  barcodes <- barcodes[sorted,]
  colsum <- sum(sorted)
  colsum.90 <- sum(sorted)*0.9
  colsum.99 <- sum(sorted)*0.99
  cumsum <- cumsum(sorted)
  cumsum <- cumsum/max(cumsum)
  length(cumsum)
  b <- seq(1,length(cumsum),1)
  b <- b/length(cumsum)
  d <- data.frame(barcode=b, proportion=cumsum)
  dim(d)

  p <- ggplot2::ggplot(d, ggplot2::aes(y = d$proportion, x = d$barcode)) +
    ggplot2::geom_point(stat = "identity", show.legend = F) +
    ggplot2::theme_bw() +
    ggplot2::scale_size_manual(values=c(2,2)) +
    ggplot2::xlab("Barcode rank (Proportion)") +
    ggplot2::ylab("Cumulative Sum (Proportion)") +
    ggplot2::ggtitle(paste(samplename, ": Barcode cumulative sum"))

  return(p)
}
