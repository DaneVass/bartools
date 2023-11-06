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
#' @param dgeObject DGEList object with barcode counts.
#' @param div Optional, precomputed diversity metrics calculated by `calcDivIndexes()`.
#' @param group Optional, column name in sample metadata to group samples by (string).
#' @param metric Diversity metric to plot (string). One of "shannon", "simpson", "invsimpson" or "gini". Default = `shannon`.
#' @param type Plot as bar graph or point (string). Default = `bar`.
#'
#' @return Returns a plot of calculated diversity index per sample
#' @importFrom vegan diversity
#' @importFrom ineq Gini
#' @export
#' @examples
#' plotDivIndexes(counts = test.dge$counts)

plotDivIndexes <-
  function(dgeObject,
           div = NULL,
           metric = "shannon",
           type = "bar",
           group = NULL) {
    inputChecks(dgeObject, groups = group)

    # check metric
    if (!metric %in% c("shannon", "simpson", "invsimpson", "gini")) {
      stop("Metric argument must be one of shannon, simpson, invsimpson, or gini.")
    }
    if (!type %in% c("bar", "point")) {
      stop("Type must be one of bar or point.")
    }

    if (is.null(div)) {
      # calculate div indices
      div <- calcDivIndexes(dgeObject)
    }

    if (!is.null(group)) {
      div[group] <- dgeObject$samples[[group]]
    }

    # barplot
    if (type == "bar") {
      p <-
        ggplot2::ggplot(div, aes_string(x = "name", y = metric, fill = group)) +
        ggplot2::geom_bar(stat = "identity")
    }

    if (type == "point") {
      p <-
        ggplot2::ggplot(div, aes_string(x = "name", y = metric, color = group)) +
        ggplot2::geom_point(stat = "identity", size = 3)
    }

    p <- p +
      ggplot2::scale_color_manual(values = rev(ggpubr::get_palette("npg", length(
        unique(div$group)
      )))) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour =
                                                                  "grey70")) +
      ggplot2::labs(title = (paste(
        "Diversity Index:", as.character(metric)
      ))) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab(as.character(metric)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (is.null(group)) {
      p <- p +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(NULL)
    }

    return(p)
  }
