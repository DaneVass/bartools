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
#' @param N String of sample 1 name (used to ranked barcode abundance)
#' @param name String of the heatmap scale name
#' @param show_bc Boolean to show barcode names on the heatmap
#'
#' @return Returns a heatmap
#' @export
#' @examples
#' data(test.counts)
#' plotBarcodeHeatmap(test.counts, 10,'Counts',FALSE)

plotBarcodeHeatmap <- function(counts, N, name = "CPM", show_bc = FALSE){
    bc <- c()
    for (s in colnames(counts)){
        top_bc = utils::head(rownames(counts[order(counts[,s],decreasing = T),]), N)
        bc <- c(bc, top_bc)
    }
    tab = counts[rownames(counts) %in% bc,]
    suppressMessages(ComplexHeatmap::Heatmap(as.matrix(tab),
                                             name = name,show_row_names = show_bc,
                                             cell_fun = function(j, i, x, y, w, h, fill) {
                                                 if(tab[i, j] %in% utils::head(sort(tab[,j], decreasing=TRUE), n=N)) {
                                                     grid::grid.text("*", x, y, vjust = 0.8)
                                                 }}))
}
