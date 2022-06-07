#' Plot sample read counts
#'
#' Simple plot of total read counts per sample
#'
#' @param counts data.frame of barcode count x sample
#' @param group a character vector of containing grouping information. Must be equal the number of columns of counts. Can pass a metadata column from DGEList object.
#' @param log10 Boolean. log10 transform output?
#' @param legend Boolean. Include legend?
#' @param order Boolean. Order samples by group?
#'
#' @return Returns a plot of the read counts per column (sample) in a data frame
#' @export
#'
#' @examples
#' data(test.counts)
#' plotReadCounts(test.dge$counts, group = test.dge$samples$Group)
#'

plotReadCounts <- function(counts, group = NULL, log10 = FALSE, legend = TRUE, order = TRUE){

  if(!is.null(group)){
    cols <- as.factor(group)
  } else {
    cols <- rep("Group 1", length(colnames(counts)))
  }

  dat <- colSums(counts)
  dat <- data.frame(sample = names(dat), counts = dat, group = cols, row.names = NULL)

  if(order){
    dat$sample <- factor(dat$sample, levels = dat$sample[order(dat$group, decreasing = T)])
  }

  # plot data
  p <- ggplot2::ggplot(data=dat, ggplot2::aes(y=sample, x=counts, fill = group)) +
    ggplot2::geom_bar(stat="identity", width = .75) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0, vjust=0, hjust = 0.5)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5)) +
    ggplot2::ggtitle("Total Read Counts") +
    ggplot2::xlim(0,max(dat$counts+100)) +
    ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(unique(dat$group)))))

  if(isTRUE(log10)){
    p <- p + ggplot2::scale_x_log10()
  }

  if(legend == FALSE){
    p <- p + ggplot2::theme(legend.position = "none")
  }

  print(p)
}
