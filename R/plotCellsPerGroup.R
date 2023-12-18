#' plotCellsPerGroup
#'
#' Plots number of cells for each detected barcode in a single cell dataset
#'
#' @title
#' Plot of cells per barcode
#'
#' @description
#' Takes a single cell object, a grouping variable, a factor within the group to test,
#' and an ident class (i.e. clusters). Per level of ident performs hypergeometric testing
#' for enrichment of group factor
#'
#' @param sc.obj single cell object in Seurat or SingleCellExperiment format containing group metadata
#' @param group a column of metadata in sc.obj
#' @param order Logical. Rank order levels of group by cell number
#' @param trans From ggplot2. For continuous scales, the name of a transformation object or the object itself.
#' @param threshold threshold number of cells above which labels for group levels will appear
#' @param plot Logical. Plot results or return data.
#' @param label Logical. Label group levels above cell number threshold. if ggrepel is installed will use geom_text_repel
#' instead of geom_text
#' @param sep Separating character used for aggregation (string). Default = `;`.
#'
#' @return Returns a plot of cell number by barcode test results or underlying plot data
#' @importFrom rlang .data
#' @export
#'
#'

plotCellsPerGroup <- function(sc.obj = NULL,
                              group = NULL,
                              order = TRUE,
                              trans = NULL,
                              threshold = 100,
                              plot = TRUE,
                              label = TRUE,
                              sep = ";") {
  # check inputs
  if (is.null(sc.obj)) {
    stop("Please supply a Seurat or SingleCellExperiment object")
  }

  # get metadata and ident class
  if (class(sc.obj)[1] == "Seurat") {
    meta <- sc.obj@meta.data
    type <- "Seurat"

  } else {
    if (class(sc.obj)[1] == "SingleCellExperiment") {
      meta <- sc.obj@colData
      type <- "SingleCellExperiment"

    } else {
      stop("Single cell object must be supplied in Seurat or SingleCellExperiment format")
    }
  }

  # extract group variable and remove NA rows
  meta.bc <- meta[!is.na(meta[, `group`]), ]

  # tally instances of each level of group
  bc.tally <- as.data.frame(table(meta.bc[, `group`]))

  if (order) {
    bc.tally <- bc.tally[order(bc.tally$Freq, decreasing = F), ]
    bc.tally$Var1 <-
      factor(bc.tally$Var1, levels = bc.tally$Var1[order(bc.tally$Freq, decreasing = T)])
  }

  # count number of barcodes
  bc.tally$num.barcodes <-
    factor(unlist(lapply(
      strsplit(
        as.character(bc.tally$Var1),
        split = sep,
        perl = T
      ), length
    )))

  # add labelling variable according to defined cell number threshold
  bc.tally$label <-
    ifelse(bc.tally$Freq > threshold, as.character(bc.tally$Var1), "")


  # plot p value histogram or return raw data
  if (plot) {
    p <- ggplot2::ggplot(bc.tally) +
      ggplot2::geom_point(ggplot2::aes(x = .data$Var1, y = .data$Freq, color = .data$num.barcodes)) +
      ggplot2::xlab("barcode") +
      ggplot2::theme_classic() +
      ggplot2::ylab("number of cells") +
      ggplot2::geom_hline(yintercept = threshold, color = "red") +
      ggplot2::scale_color_manual(values = viridis::viridis(length(unique(
        bc.tally$num.barcodes
      )))) +
      ggplot2::ggtitle("Number of cells per clone")

    if (!is.null(trans)) {
      p <- p + ggplot2::scale_y_continuous(trans = trans)
    }

    if (label) {
      if (rlang::is_installed("ggrepel")) {
        p <- p + ggrepel::geom_text_repel(
          data = bc.tally,
          ggplot2::aes(
            x = .data$Var1,
            y = .data$Freq,
            color = .data$num.barcodes,
            label = .data$label
          ),
          max.overlaps = 20
        )

      } else {
        p <- p + ggplot2::geom_text(
          data = bc.tally,
          ggplot2::aes(
            x = .data$Var1,
            y = .data$Freq,
            color = .data$num.barcodes,
            label = .data$label
          ),
          nudge_x = -20
        )

      }
    }
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
    return(p)

  } else {
    colnames(bc.tally) <- c("Barcode", "Freq", "Num.Barcodes", "Label")
    return(bc.tally)
  }
}
