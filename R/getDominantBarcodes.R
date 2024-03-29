#' @title
#' Get dominant barcodes per sample
#'
#' @description
#' Takes a DGEList or dataframe of barcode counts,
#' computes the proportion abundance of each barcode within each sample and
#' then returns a list of barcodes meeting a threshold abundance per sample.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param threshold Proportion threshold to call dominant barcodes (decimal). Default = `0.05` (i.e. 5 percent).
#'
#' @return Returns a named list containing vectors of dominant barcodes per sample
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' data(test.dge)
#' getDominantBarcodes(test.dge, threshold = 0.05)

getDominantBarcodes <- function(dgeObject, threshold = 0.05) {
  inputChecks(dgeObject)

  # transform counts into proportion within sample
  barcodes.proportional <- as.data.frame(dgeObject$counts)
  barcodes.proportional <-
    sweep(barcodes.proportional,
          2,
          colSums(barcodes.proportional),
          `/`)

  # select all barcodes above threshold per sample
  barcodes.proportional.dominant <- barcodes.proportional %>%
    tibble::rownames_to_column("barcode") %>%
    tidyr::pivot_longer(-.data$barcode, names_to = "sample", values_to = "proportion") %>%
    dplyr::filter(.data$proportion > threshold) %>%
    dplyr::arrange(dplyr::desc(.data$proportion))

  samps <- colnames(barcodes.proportional)

  dominant.barcodes <- lapply(samps, function(x) {
    barcodes.proportional.dominant %>%
      dplyr::filter(.data$sample == x) %>%
      dplyr::pull(.data$barcode)
  })
  names(dominant.barcodes) <- samps

  return(dominant.barcodes)
}
