#' plotMetrics
#'
#' Plots metric distribution for selected groups in a single cell dataset
#'
#' @title
#' Violin plot of selected metric per group variable
#'
#' @description
#' Takes a single cell object, a grouping variable and a metadata variable.
#' Plots distribution of factor per level of group with total number of cells above threshold
#'
#' @param sc.obj single cell object in Seurat or SingleCellExperiment format containing group metadata
#' @param group a column of metadata in sc.obj containing grouping information
#' @param trans From ggplot2. For continuous scales, the name of a transformation object or the object itself.
#' @param factor a metadata field to plot per level of group
#' @param threshold threshold number of cells per level of group to plot
#' @param plot Logical. Plot results or return data. 
#'
#' @return Returns a plot of cell number by barcode test results or underlying plot data
#' @export
#' @examples
#' 
#' 

plotMetrics <- function(sc.obj = NULL, 
                                group = NULL,
                                trans = NULL,
                                factor = NULL,
                                threshold = 100,
                                plot = TRUE) {
  
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
      stop("A single cell object must be supplied in Seurat or SingleCellExperiment format")
    }  
  }
  
  # extract group metadata
  meta.bc <- as.data.frame(meta[!is.na(meta[,`group`]),])
  bc.tally <- as.data.frame(table(meta.bc[,`group`]))
  select <- (bc.tally[which(bc.tally$Freq > threshold),])$Var1
  meta.select <- meta.bc[which(meta.bc[,`group`] %in% select),]
  meta.select.group <- dplyr::group_by(meta.select, barcode)
  
  # convert grouping var to factor if not already
  meta.select[,`group`] <- factor(meta.select[,`group`])
  
  # summarise per clone
  clone.summary <- dplyr::summarise_at(meta.select.group, vars(factor), 
                                       list(avg = mean, max = max, 
                                            min = min, sd = sd))
  
  # plot distributions (optional box or violin plot) 
  if (plot) {
    group <- ggplot2::sym(group)
    factor <- ggplot2::sym(factor)
    p <- ggplot2::ggplot(meta.select, 
                         ggplot2::aes(y = !!group, x = !!factor)) +
      ggplot2::geom_jitter(width = 0.1, size = 0.2) +
      ggplot2::geom_violin(alpha = 0.8) + 
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste(as.character(factor), "per barcode"))
    
    if (!is.null(trans)) {
      p <- p + ggplot2::scale_x_continuous(trans = trans)
    }
    return(p)
  } else {
    return(meta.select)
  }
}