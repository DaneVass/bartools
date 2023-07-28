## plotCellsInClusters
#'
#' Plot number or percentage of a group per level of ident in a single
#' cell object
#'
#' @title
#' Plot cells per ident class
#'
#' @description
#' Takes a single cell object, a grouping variable, a factor within the group to test, 
#' and an ident class (i.e. clusters). Plots percentage or raw number of cells within factor
#' present per level of ident (i.e. per cluster)
#'
#' @param sc.obj single cell object in Seurat or SingleCellExperiment format containing idents and group metadata
#' @param group a column of metadata 
#' @param factor a level of group to plot per level of idents
#' @param idents a column of metadata defining the cluster identities of each cell. default = "seurat_clusters"
#' @param plot.pct Logical. Plot percentages (TRUE) or raw numbers (FALSE). default = TRUE
#' @param plot Logical. Plot or return raw data. default = TRUE
#'
#' @return Returns a histogram or underlying plot data
#' @export


plotCellsInClusters <- function(sc.obj,
                                group = NULL, 
                                factor = NULL,
                                idents = NULL, 
                                plot.pct = T,
                                plot = T){
  
  # check inputs
  if (is.null(sc.obj)) {
    stop("Please supply a Seurat or SingleCellExperiment object")
  }
  
  # get metadata and ident class
  if (class(sc.obj)[1] == "Seurat") {
    meta <- sc.obj@meta.data
    type <- "Seurat"
    if (is.null(idents)) {
      idents <- "seurat_clusters"
    }
  } else {
    if (class(sc.obj)[1] == "SingleCellExperiment") {
      meta <- sc.obj@colData
      type <- "SingleCellExperiment"
      if (is.null(idents)) {
        idents <- "ident"
      }
    } else {
      stop("Single cell object must be supplied in Seurat or SingleCellExperiment format")
    }  
  }
  
  # get number of seurat clusters. NB plotting is 0 indexed so total - 1
  clusters <- data.frame(clusters = as.factor(seq(0,max(as.numeric(meta[,idents])-1), 1)))
  colnames(clusters) <- as.character(idents)
  # get and count cells of interest
  meta.bc <- as.data.frame(meta[!is.na(meta[,`group`]),])
  # get cells in group that are factor
  meta.select <- meta.bc[which(meta.bc[,`group`] == factor),]
  
  # group by and tally clusters (i,e. idents)
  idents <- sym(idents)
  meta.bc.group <- dplyr::tally(dplyr::group_by(meta.select, !!idents))
  print(meta.bc.group)
  
  # join counts per cluster (ident) of group factor into table of clusters (idents)
  dat <- dplyr::left_join(clusters, meta.bc.group, by = "seurat_clusters")
  dat[is.na(dat)] = 0
  
  # convert to proportion
  dat <- dplyr::mutate(dat, Percentage = 100*(n/sum(dat$n)))
  
  # plot data
  if(plot){
    if(isTRUE(plot.pct)){
      p <- ggplot(dat, aes(x = !!idents, y = Percentage, fill = !!idents)) +
        geom_histogram(stat = "identity") +
        scale_fill_manual(values=scales::hue_pal()(nrow(clusters))) +
        theme_bw() +
        ggtitle(paste("Percentage of cells in clusters:", group, factor))
    } else {
      p <- ggplot(dat, aes(x = !!idents, y = n, fill = !!idents)) +
        geom_histogram(stat = "identity") +
        scale_fill_manual(values=scales::hue_pal()(nrow(clusters))) +
        theme_bw() +
        ggtitle(paste("Number of cells in clusters:", group, factor))
    }  
    return(p)
  } else {
    return(dat)
  }
}