#' Annotate Barcode Reads
#'
#' Matches cells in a single cell experiment to detected DNA barcodes.
#' Dataframe returned will have all cells matched to a barcode
#' If there is no barcode matchable to a cell "not.detected" is returned
#' For cells that have multiple detected barcodes each barcode is returned separated by ';'
#'
#' @param all.cells list of all cells in the experiment to be annotated
#' @param barcoded.cells dataframe a dataframe with cell id and barcode id as columns
#'
#'
#' @return Returns a data-frame containing the the 10X cell ID and the lintraceR DNA barcode ID.
#'
#' @import magrittr
#' @export
#'

annotateBarcodeReads <- function(all.cells, barcoded.cells){

  message("Annotating 10X cell ID with DNA barcode ID")
  # make sure this takes the cell id and barcode reference id columns
  cells.w.barcode <- dplyr::select(barcoded.cells$Cell.10X.Barcode, barcoded.cells$referenceID, .data = barcoded.cells)
  anno.merge <- dplyr::left_join(all.cells, cells.w.barcode, by=c("V1"="Cell.10X.Barcode"))

  # order by cell id and barcode number within cell and get unique entries
  vals <- as.numeric(gsub("[A-Z]*_Barcode_","", anno.merge$referenceID, ignore.case = T))
  anno.order <- anno.merge[order(anno.merge$V1,vals),]
  anno.uniq <- anno.order[!duplicated(anno.order),]

  # collapse multiple barcodes per cell
  anno.collapse <- stats::aggregate(anno.uniq$referenceID, list(anno.uniq$V1), paste, collapse=";")

  # reformat
  anno.collapse$x <- gsub("NA", "not.detected", anno.collapse$x)
  colnames(anno.collapse) <- c("cell.id", "barcode")
  anno.final <- data.frame(barcode = anno.collapse$barcode, row.names = gsub('-1', '', anno.collapse$cell.id))

  # Annotation summary stats
  print("Total cells annotated:")
  print(nrow(anno.final))

  print("Total cells with barcode:")
  print(nrow(anno.final %>% dplyr::filter(anno.final$barcode != "not.detected")))

  print("Total cells without barcode:")
  print(nrow(anno.final %>% dplyr::filter(anno.final$barcode == "not.detected")))

  print("Percentage of cells with annotated barcode:")
  print(nrow(anno.final %>% dplyr::filter(anno.final$barcode != "not.detected"))/nrow(anno.final)*100)

  print("Number of unique barcodes:")
  print(length(unique(anno.final$barcode)))

  return(anno.final)
}

