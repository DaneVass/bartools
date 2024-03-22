#' @title
#' aggregateBarcodesClones
#'
#' @description
#' Aggregate barcodes and barcode UMI counts into clones and subsequently aggregate clones per cell, all in alphabetical order.
#'
#' @param counts Count table with columns cellid, barcode and bc.umi.count (output of `readBartabCounts()`).
#' @param clones A named list with barcodes as names and clone identity as values (output of `barcodesToClones()`).
#' @param sep Separator between barcodes (string). Default = `;`.
#'
#' @return Returns a data frame with one row per cell ID, clone ID, barcode ID and barcode UMI counts.
#'
#' @importFrom rlang .data
#' @export
#

aggregateBarcodesClones <-
  function(counts,
           clones,
           sep = ";") {

    counts <- m6_clone_counts_mch
    clones <- m6_clones
    counts <- counts %>%
      dplyr::mutate(clone_id = paste0("clone_", clones[barcode]))

    # make sure barcodes are aggregated in alphabetical order to ensure comparability
    counts <- dplyr::arrange(counts, .data[["clone_id", "barcode"]])
    bc.counts <- as.data.frame(
      data.table::dcast(
        data.table::setDT(counts),
        cellid ~ .,
        value.var = c("clone_id", "barcode", "bc.umi.count"),
        fun.aggregate = function(x) {
          # unique clones but not counts
          if (is.numeric(x))
            paste(x, collapse = sep)
          else
            paste(unique(x), collapse = sep)
        }
      )
    )
    return(bc.counts)
  }
