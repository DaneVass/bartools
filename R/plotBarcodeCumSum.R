#' plotBarcodeCumSum
#'
#' Plot the cumulative sum of barcode cpm counts for a list of samples against and ordered by sample1.
#'
#' @title
#' Barcode cumulative sum plot of a list of samples against and ordered by sample1
#'
#' @description
#' Takes a dataframe of the barcode cpm counts, calculate the relative abundance of each barcode in each sample.
#' Then, barcodes of all samples are ordered for sample1 in decreasing order and the cumulative sum is calculated for
#' each sample. Cumulative sum is then plotted against sample1.
#'
#' @param counts Dataframe containing cpm counts
#' @param sample1 sample to compare others against
#' @param samples vector of sample names to be plotted against sample1
#'
#' @return Returns a plot of the cumulative sum of read counts per barcode
#' @export
#'
#' @examples
#' data(test.counts)
#' plotBarcodeCumSum(counts = test.counts, sample1 = 'T0-1', samples = c('T0-1','T0b-1'))

plotBarcodeCumSum <- function(counts, sample1 = NULL, samples = NULL) {

  var1 <- counts[,sample1]
  abundance1 <- var1/sum(var1)

  plot(1,
       xlab = paste('Cumulative sum abundance of ',sample1,sep=''), ylab= 'Cumulative sum abundance',
       xlim = c(0,1), ylim = c(0,1),
       main= paste('Cumulative abundance ranked by ',sample1,sep='') )
  cols <- ggpubr::get_palette("npg", length(samples))
  i = 1
  for (s in samples){
    var2 <- counts[,s]
    abundance2 <- var2/sum(var2)
    df <- data.frame(abundance1,abundance2)
    df <- df[order(var1,decreasing = T),]
    df$abundance1 <- cumsum(df$abundance1)
    df$abundance2 <- cumsum(df$abundance2)
    colnames(df) <- c('CumSum1','CumSum2')
    # add 0's in the first for plotting
    x <- rep(0, ncol(df))
    df <- rbind(x, df)
    graphics::lines(df$CumSum1,df$CumSum2,type='l',lwd=3,col = cols[i])
    i=i+1
  }
  graphics::legend(0,1,legend=samples,col = cols,lty=1)
  graphics::abline(a=0,b=1,col='black')
}

#
#
#
#   # filter requested samples
#   if(!is.null(samples)){
#     counts <- counts[,which(samples %in% colnames(counts)), drop = F]
#     colnames(counts) == samples
#   }
#
#   # setup output dataframe
#   cumsum.df <- data.frame(row.names = seq(1:nrow(counts)))
#
#   # sort barcodes in samples
#   for (samp in colnames(counts)){
#     ordered <- order(counts[,samp], decreasing = T)
#     sorted <- counts[ordered,samp,drop = F]
#     colsum <- sum(sorted)
#     percentile <- sum(sorted)*pct
#     cumsum <- cumsum(sorted)
#     cumsum <- cumsum/max(cumsum)
#
#     b <- seq(1,length(cumsum),1)
#     b <- b/length(cumsum)
#     d <- data.frame(proportion=cumsum)
#     dim(d)
#
#     cumsum.df <- cbind(cumsum.df, d)
#   }
#
#   # plot cumsums
#   cumsum.df <- reshape2::melt(cumsum.df)
#   ggplot2::ggplot(cumsum.df, ggplot2::aes(1:nrow(counts), y=value, color = variable)) +
#     ggplot2::geom_point(stat = "identity", show.legend = F) +
#     ggplot2::geom_smooth(stat = "identity", show.legend = F) +
#     ggplot2::theme_bw() +
#     ggplot2::scale_size_manual(values=c(2,2)) +
#     ggplot2::xlab("Barcode rank (Proportion)") +
#     ggplot2::ylab("Cumulative Sum (Proportion)") +
#     ggplot2::ggtitle("Barcode cumulative sum")

