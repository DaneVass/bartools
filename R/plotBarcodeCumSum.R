#' plotBarcodeCumSum
#'
#' Plot the cumulative sum of barcode CPM counts for a list of samples against and ordered by sample1.
#'
#' @title
#' Barcode cumulative sum plot of a list of samples against and ordered by reference sample
#'
#' @description
#' Takes a dataframe of the barcode cpm counts, calculate the relative abundance of each barcode in each sample.
#' Then, barcodes of all samples are ordered for the reference sample in decreasing order and the cumulative sum is calculated for
#' each sample. Cumulative sum is then plotted against the reference sample
#'
#' @param dgeObject DGE object containing CPM counts
#' @param referenceSample sample to compare others against
#' @param samples vector of sample names to be plotted against reference sample
#'
#' @return Returns a plot of the cumulative sum of read counts per barcode
#' @export
#'
#' @examples
#' data(test.dge)
#' plotBarcodeCumSum(test.dge, referenceSample = 'T0-1', samples = c('T0-1','S9-1', 'S8-1'))
plotBarcodeCumSum <-
  function(dgeObject,
           referenceSample = NULL,
           samples = NULL) {
    inputChecks(dgeObject, samples = c(samples, referenceSample))
    counts <- as.data.frame(dgeObject$counts)

    var1 <- counts[, referenceSample]
    abundance1 <- var1 / sum(var1)

    plot(
      1,
      xlab = paste('Cumulative sum abundance of ',
                   referenceSample, sep =
                     ''),
      ylab = 'Cumulative sum abundance',
      xlim = c(0, 1),
      ylim = c(0, 1),
      main = paste('Cumulative abundance ranked by ', referenceSample, sep =
                     '')
    )
    cols <- colors.npg[1:length(samples)]
    i = 1
    for (s in samples) {
      var2 <- counts[, s]
      abundance2 <- var2 / sum(var2)
      df <- data.frame(abundance1, abundance2)
      df <- df[order(var1, decreasing = T),]
      df$abundance1 <- cumsum(df$abundance1)
      df$abundance2 <- cumsum(df$abundance2)
      colnames(df) <- c('CumSum1', 'CumSum2')
      # add 0's in the first for plotting
      x <- rep(0, ncol(df))
      df <- rbind(x, df)
      graphics::lines(
        df$CumSum1,
        df$CumSum2,
        type = 'l',
        lwd = 3,
        col = cols[i]
      )
      i = i + 1
    }
    graphics::legend(0,
                     1,
                     legend = samples,
                     col = cols,
                     lty = 1)
    graphics::abline(a = 0, b = 1, col = 'black')
  }
