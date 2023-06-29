#' plotBCBoxplot
#'
#' Boxplot of selected barcodes
#'
#' @title
#' Boxplot of selected barcodes
#'
#' @description
#' Takes a dataframe of barcode counts,
#' select specific barcodes in the dataframe
#' and boxplot their counts
#'
#' @param dge.obj DGE object containing counts and sample information
#' @param barcodes a vector of barcodes to be plotted
#' @param group metadata field to facet data on
#' @param conditions specific levels of group to be plotted. If NULL all levels are plotted
#' @param trans From ggplot2. For continuous scales, the name of a transformation object or the object itself.
#' @param point From ggplot2. Include plot points?
#' @param violin From ggplot2. Include violin plot?
#'
#' @return Returns a boxplot
#' @export
#' @examples
#' data(test.dge)
#' plotBarcodeBoxplot(test.dge, barcodes = "BC_190202", group = "Treatment", 
#' conditions = c("Vehicle", "Low_dose", "High_dose"))

plotBarcodeBoxplot <- function(dge.obj, 
                               barcodes = NULL,
                               group = NULL,
                               conditions = NULL,
                               trans = NULL, 
                               point = FALSE,
                               violin = FALSE){
  
  # check inputs
  if (methods::is(dge.obj)[1] == "DGEList") {
    counts.raw <- as.data.frame(dge.obj$counts)
  }
  else {
    stop("input must be a DGEList object")
  }
  
  if(is.null(barcodes)){
    stop("please provide some barcodes to plot")
  }
    
  sample = dge.obj$samples
  counts.cpm <- cpm(counts.raw)
  
  # no metadata group specified
  if(is.null(group)){
    
    # only one barcode given
    if (length(barcodes) == 1){
      message("single barcode given. No metadata group")
      dat <- counts.cpm[which(rownames(counts.cpm) == barcodes),]
      p <- ggplot2::ggplot() + 
        ggplot2::geom_boxplot(aes(y = dat)) + 
        ggplot2::labs(xlab = as.factor(barcodes)) +
        ggplot2::theme_bw()
      
      # check violin
      if (violin) {
        p <- ggplot2::ggplot() + 
          ggplot2::geom_violin(aes(x = barcodes, y = dat), scale = "width") +
          ggplot2::geom_boxplot(aes(x = barcodes, y = dat), width = 0.2) + 
          ggplot2::labs(xlab = as.factor(barcodes))
      }
      
      if (!is.null(trans)) { p <- p + ggplot2::scale_y_continuous(trans = trans) } # check trans
      if (point) { p <- p + ggplot2::geom_point(aes(x = barcodes, y = dat)) } # check point
      
      return(p)
      
    } 
    else {
      dat <- counts.cpm[which(rownames(counts.cpm) %in% barcodes),,drop = F]
      dat <- reshape2::melt(dat)
      colnames(dat) <- c("barcodes", "sample", "value")
      message("multiple barcodes given. No metadata group")
      
      p <- ggplot2::ggplot(dat) + 
        ggplot2::geom_boxplot(aes(x = barcodes, y = value, fill = barcodes)) + 
        ggplot2::xlab(NULL) +
        ggplot2::ylab("CPM") +
        ggplot2::theme_bw()
      
      # check violin
      if (violin) {
        p <- ggplot2::ggplot(dat) + 
          ggplot2::geom_violin(aes(x = barcodes, y = value, fill = barcodes), scale = "width") +
          ggplot2::geom_boxplot(aes(x = barcodes, y = value), width = 0.2) + 
          ggplot2::labs(xlab = as.factor(barcodes))
      }
      
      if (!is.null(trans)) { p <- p + ggplot2::scale_y_continuous(trans = trans) } # check trans
      if (point) { p <- p + ggplot2::geom_point(aes(x = barcodes, y = value)) } # check point
      
      return(p)
      
    }
    
    # plot with facets for metadata group and condition
    
  } else {
    # get group metadata
    if (length(group) > 1) { stop("please enter a single grouping variable") }
    group.col <- which(colnames(dge.obj$samples) == group)
    group.dat <- dge.obj$samples[,group.col,drop = F]
    group.dat$sample <- rownames(group.dat)
    
    # only one barcode given with a metadata group
    if (length(barcodes) == 1){
      message(paste("single barcode given. Metadata group = ", group))
      dat <- as.data.frame(counts.cpm[which(rownames(counts.cpm) == barcodes),])
      dat$sample <- rownames(dat)
      
      # incorporate metadata
      dat <- dplyr::left_join(dat, group.dat, by = "sample")
      colnames(dat) <- c(barcodes, "sample", group)
      dat[,2:length(colnames(dat))] <- lapply(dat[,2:length(colnames(dat))], factor)
      
      # subset on selected conditions
      if (!is.null(conditions)) {
        dat <- dat[which(dat[,3] %in% conditions),]
        if(nrow(dat) == 0){
          stop(paste("requested conditions", paste(conditions, collapse = " "),"were not found or were entered incorrectly. stopping"))
        }
      }
      
      p <- ggplot2::ggplot() + 
        ggplot2::geom_boxplot(data = dat, mapping = aes(y = dat[,barcodes], x = dat[,group], fill = dat[,group])) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("CPM") +
        ggplot2::theme_bw()
      
      # check violin
      if (violin) {
        p <- ggplot2::ggplot() + 
          ggplot2::geom_violin(aes(y = dat[,barcodes], x = dat[,group], fill = dat[,group]), scale = "width") +
          ggplot2::geom_boxplot(aes(y = dat[,barcodes], x = dat[,group]), width = 0.2) + 
          ggplot2::xlab(NULL) +
          ggplot2::ylab("CPM") +
          ggplot2::theme_bw()
      }
      
      if (!is.null(trans)) { p <- p + ggplot2::scale_y_continuous(trans = trans) +
        ggplot2::ylab(paste(trans, "CPM")) } # check trans
      
      if (point) { p <- p + 
        ggplot2::geom_point(aes(y = dat[,barcodes], x = dat[,group], fill = dat[,group])) } # check point
      
      return(p)
    } else {
      # multiple barcodes given (facet on the barcodes)
      message(paste("Multiple barcodes given. Metadata group = ", group))
      dat <- counts.cpm[which(rownames(counts.cpm) %in% barcodes),,drop = F]
      dat <- reshape2::melt(dat)
      colnames(dat) <- c("barcodes", "sample", "value")
      
      # incorporate metadata
      dat <- dplyr::left_join(dat, group.dat, by = "sample")
      dat[,c(1,2,4)] <- lapply(dat[,c(1,2,4)], factor)
      
      # subset on selected conditions
      if (!is.null(conditions)) {
        dat <- dat[which(dat[,4] %in% conditions),]
        if(nrow(dat) == 0){
          stop(paste("requested conditions", paste(conditions, collapse = " "),"were not found or were entered incorrectly. stopping"))
        }
      }
      
      p <- ggplot2::ggplot(dat) + 
        ggplot2::geom_boxplot(aes(y = value, x = dat[,group], fill = dat[,group])) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("CPM") +
        ggplot2::theme_bw()
        
      
      # check violin
      if (violin) {
        p <- ggplot2::ggplot(dat) + 
          ggplot2::geom_violin(aes(y = value, x = dat[,group], fill = dat[,group]), scale = "width") +
          ggplot2::geom_boxplot(aes(y = value, x = dat[,group]), width = 0.2) + 
          ggplot2::xlab(NULL) +
          ggplot2::ylab("CPM") +
          ggplot2::theme_bw()
      }
      
      if (!is.null(trans)) { p <- p + ggplot2::scale_y_continuous(trans = trans) +
        ggplot2::ylab(paste(trans, "CPM")) } # check trans
      
      if (point) { p <- p + ggplot2::geom_point(aes(y = value, x = dat[,group])) } # check point
      
      p <- p + ggplot2::labs(fill=as.character(group)) 
      
      return(p + facet_wrap(~barcodes))
      
    }
  }
}
