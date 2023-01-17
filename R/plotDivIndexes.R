#' plotDivIndexes
#'
#' Plots common diversity indices per sample
#'
#' @title
#' Plot diversity indices across barcode samples
#'
#' @description
#' Plots a dataframe of diversity indices per sample
#' 
#' @param counts dataframe of raw or normalised counts with samples as columns and observations/barcodes as rows
#' @param div dataframe of diversity metrics as columns and samples as rows
#' @param metric diversity metric to plot. one of "shannon", "simpson", "invsimpson" or "gini"
#'
#' @return Returns a plot of calculated diversity index per sample
#' @importFrom vegan diversity
#' @importFrom ineq Gini
#' @export
#' @examples
#' plotDivIndexes(dat = test.counts)

plotDivIndexes <- function(counts, div = NULL, metric = "shannon", type = "bar"){
  # calculate div indices if counts given
  if (is.null(div)){
    div <- calcDivIndexes(counts)  
  }
  
  # check metric
  if(!metric %in% c("shannon", "simpson", "invsimpson", "gini")){
    stop("metric argument must be one of shannon, simpson, invsimpson, or gini")
  }
  
  # initialise metric string for plotting
  metric <- metric
  
  # barplot
  if (type == "bar"){
    p <- ggplot2::ggplot(div, aes(x = name, y = .data[[metric]], fill = group)) +
      ggplot2::geom_bar(stat = "identity") + 
      ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(unique(div$group))))) +
      ggplot2::theme(panel.grid.major.x=ggplot2::element_line(colour="grey70")) +
      ggplot2::labs(title = (paste("Diversity Index:",as.character(metric)))) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab(as.character(metric)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    p
    
    return(p)  
  }
  
  if (type == "point"){
    p <- ggplot2::ggplot(div, aes(x = group, y = .data[[metric]], color = group)) +
      ggplot2::geom_point(stat = "identity", size = 3) + 
      ggplot2::scale_color_manual(values = rev(ggpubr::get_palette("npg", length(unique(div$group))))) +
      ggplot2::theme(panel.grid.major.x=ggplot2::element_line(colour="grey70")) +
      ggplot2::labs(title = (paste("Diversity Index:",as.character(metric)))) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab(as.character(metric)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    p
    return(p)
  }
  
}
