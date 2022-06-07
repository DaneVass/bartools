#' Threshold counts
#'
#' Threshold dataframe to a given level and return number of barcodes meeting threshold in each
#' sample.
#'
#' @param df Dataframe to be thresholded.
#' @param threshold The threshold to use. Rows with count below this will be removed.
#' @param plot Logical. Draw plots of dataset?
#' @return Returns a thresholded data-frame
#'
#' @export
#' @examples
#' thresholdCounts(test.counts, threshold = 20, plot = FALSE)

thresholdCounts <- function(df, threshold = 20, plot = FALSE){

  if(!is.null(threshold)){
    threshold <- threshold
  } else {
    message("No threshold given, defaulting to 20.")
    threshold <- 20
  }

  # setup output dataframe
  above.threshold.counts <- data.frame(Sample=factor(), Count=c())

  for(sample in colnames(df)){
    above.threshold = length(which(df[,sample] >= threshold))
    d <- data.frame(Sample=factor(sample),Count=above.threshold)
    above.threshold.counts <- rbind(above.threshold.counts, d)
  }

  if(plot == FALSE){
    return(above.threshold.counts)
  } else {
    g <- ggplot2::ggplot(above.threshold.counts, ggplot2::aes(x=above.threshold.counts$Sample,y=above.threshold.counts$Count))
    g + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme(panel.grid.major.x=ggplot2::element_line(colour="grey70")) +
      ggplot2::labs(title = paste("Number of barcodes meeting threshold:", threshold)) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("Number of barcodes") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
}
