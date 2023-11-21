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
#' @param type Plot as bar, point or box plot (string). Default = `bar`.
#'
#' @return Returns a plot of calculated diversity index per sample
#' @importFrom vegan diversity
#' @importFrom ineq Gini
#' @export
#' @examples
#' data(test.dge)
#' plotDivIndexes(counts = test.dge)

plotDivIndexes <-
  function(dgeObject,
           div = NULL,
           metric = "shannon",
           type = "bar",
           group = NULL,
           color = NULL) {
    inputChecks(dgeObject, groups = group)

    # check metric
    if (!metric %in% c("shannon", "simpson", "invsimpson", "gini")) {
      stop("Metric argument must be one of shannon, simpson, invsimpson, or gini.")
    }
    if (!type %in% c("bar", "point", "box")) {
      stop("Type must be one of bar or point.")
    }

    if (is.null(div)) {
      # calculate div indices
      div <- calcDivIndexes(dgeObject)
    }

    if (!is.null(group)) {
      div[group] <- dgeObject$samples[[group]]
    }
    if (is.null(group) & type == "box") {
      stop("Provide a group variable for a box plot.")
    }

    # barplot
    if (type == "bar") {
      p <-
        ggplot2::ggplot(div, aes_string(x = "name", y = metric, fill = group)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(
          unique(div[[group]])
        ))))
    }

    if (type == "point") {
      p <-
        ggplot2::ggplot(div, aes_string(x = "name", y = metric, color = group)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_color_manual(values = rev(ggpubr::get_palette("npg", length(
          unique(div[[group]])
        ))))
    }

    if (type == "box") {
      p <-
        ggplot2::ggplot(div, aes_string(x = group, y = metric, fill = group)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(
          unique(div[[group]])
        ))))
    }

    p <- p +
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
