#' plotBarcodeHeatmap
#'
#' Plot a heatmap of the N most abundant barcodes per sample
#'
#' @title
#' Plot barcode heatmap
#'
#' @description
#' Takes a dataframe of the barcode counts, and selects the N most abundant barcodes per sample.
#' Most abundant barcodes are gathered in a list and a heatmap for those barcodes is plotted for all samples.
#' Stars indicate the most abundant barcodes per sample.
#'
#' @param counts Dataframe containing cpm counts
#' @param N Number of top barcodes per sample
#' @param name String of the heatmap scale name
#' @param show_bc Boolean to show barcode names on the heatmap
#' @param samples Dataframe containing sample metadata sheet
#' @param group metadata fields to annotate samples
#' @param col_annot list of vectors assigning colours to each level of metadata
#' @param discrete Boolean to show presence/absence of barcodes instead of abundance
#' @param discrete_threshold Threshold for presence of barcode in samples
#'
#' @return Returns a heatmap
#' @export
#' @examples
#' plotBarcodeHeatmap(counts = cpm(test.dge$counts),N = 5,show_bc = TRUE,samples = test.dge$samples)

plotBarcodeHeatmap <- function(counts, N, name = "CPM", show_bc = FALSE, 
                               samples = NULL, group = NULL,col_annot = NULL, 
                               discrete = F,discrete_threshold = 1){
    set.seed(1)
    bc <- c()
    # select top N barcodes for each sample
    for (s in colnames(counts)){
        top_bc = utils::head(rownames(counts[order(counts[,s],decreasing = T),]), N)
        bc <- c(bc, top_bc)
    }
    tab = counts[rownames(counts) %in% bc,]
    # add annotation of samples
    if (!is.null(samples) && !is.null(group)) {
      # check that metadata names are in samples
        if (!all(group %in% colnames(samples))) {
            stop("group must be column in samples")
        }
      # keep group order for legend in plot
      group_cols <- group[group %in% colnames(samples)]
      dff <- samples[, group_cols, drop = FALSE]
      colnames(dff) <- group_cols
      
      ha <- ComplexHeatmap::HeatmapAnnotation(df = dff,
                                                  which = 'col',
                                                  col = col_annot,
                                                  annotation_width = unit(c(1, 4), 'cm'),
                                                  gap = unit(1, 'mm'))
    } else {
        ha <- NULL
    }
    if (discrete){
      tab2 <- ifelse(counts[rownames(counts) %in% bc, ] > discrete_threshold, 1, 0)
      col <- c("1" = "red", "0" = "blue")
      
      suppressMessages(ComplexHeatmap::Heatmap(
        as.matrix(tab2),
        name = "Presence",
        show_row_names = show_bc,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (tab[i, j] %in% utils::head(sort(tab[, j], decreasing = TRUE), n = N)) {
            grid::grid.text("*", x, y, vjust = 0.8)
          }
        },
        top_annotation = ha,
        col = col,
        heatmap_legend_param = list(labels = c("Yes","No"))
      ))
    } else{
        suppressMessages(ComplexHeatmap::Heatmap(as.matrix(tab),
                                             name = name,show_row_names = show_bc,
                                             cell_fun = function(j, i, x, y, w, h, fill) {
                                                 if(tab[i, j] %in% utils::head(sort(tab[,j], decreasing=TRUE), n=N)) {
                                                     grid::grid.text("*", x, y, vjust = 0.8)
                                                 }},
                                             top_annotation = ha))
    }
}
