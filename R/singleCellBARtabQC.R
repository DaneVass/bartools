#' readBartabCounts
#'
#' Read BARtab counts table from single-cell sample
#'
#' @param countsPath Path to BARtab counts.tsv file
#'
#' @return Returns a data frame in long format of cell IDs, barcodes and barcode UMI counts
#'
#' @export
readBartabCounts <- function(countsPath) {
  # read in barcode-cell pairs
  counts <- utils::read.delim(countsPath, header = TRUE)
  if (ncol(counts) != 3) {
    stop("Data is not in the expected format. Expects columns <gene, cell, count>.")
  }
  if (any(colnames(counts) != c("gene", "cell", "count"))) {
    stop("Data is not in the expected format. Expects columns <gene, cell, count>.")
  }
  counts <- counts[, c("cell", "gene", "count")]
  colnames(counts) <- c("cellid", "barcode", "bc.umi.count")
  return(counts)
}

#' aggregateBarcodes
#'
#' Aggregate barcodes and barcode UMIs per cell, in alphabetical order.
#'
#' @param counts Dataframe with barcodes and UMI counts per cell
#' @param sep Separator between barcodes (string). Default = `;`.
#'
#' @return Returns a data frame with one row per cell ID
#'
#' @importFrom rlang .data
#' @export
#
aggregateBarcodes <- function(counts, sep = ";") {
  # make sure barcodes are aggregated in alphabetical order to ensure comparability
  counts <- dplyr::arrange(counts, .data$barcode)
  bc.counts <- as.data.frame(data.table::dcast(
    data.table::setDT(counts),
    cellid ~ .,
    value.var = c("barcode", "bc.umi.count"),
    fun.aggregate = function(x) {
      paste(x, collapse = sep)
    }
  ))
  return(bc.counts)
}

#' filterBarcodes
#'
#' Filter barcodes from single-cell sample
#'
#' @param counts Dataframe with barcodes and UMI counts per cell
#' @param umiCountFilter Minimum number of UMIs per barcode per cell
#' @param umiFractionFilter Minimum fraction of UMIs per barcode per cell compared to dominant barcode in cell (barcode supported by most UMIs)
#'
#' @return Returns a data frame with one row per cell ID
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export
#

filterBarcodes <-
  function(counts,
           umiCountFilter = 1,
           umiFractionFilter = 0.3) {
    counts_filtered <- counts %>%
      dplyr::filter(.data$bc.umi.count >= umiCountFilter)

    counts_filtered <- counts_filtered %>%
      dplyr::group_by(.data$cellid) %>%
      dplyr::mutate(max_count = max(.data$bc.umi.count)) %>%
      dplyr::arrange(.data$cellid) %>%
      dplyr::filter(.data$bc.umi.count / umiFractionFilter >= .data$max_count)

    return(counts_filtered)
  }


integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

#' plotBarcodesPerCell
#'
#' Plot number of detected barcodes per cell.
#'
#' @param counts Dataframe with barcodes and UMI counts per cell or SingleCellExperiment or Seurat object
#' @param fraction Boolean, whether to plot fraction or number of cells
#' @param aggregated Counts were aggregated per cell (boolean). Default = `FALSE`.
#' @param sep Separating character used for aggregation (string). Default = `;`.
#' @param notDetected Optional, string representing no detected barcode. NA is always treated as not detected.
#'
#' @return Returns a plot
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export

plotBarcodesPerCell <- function(counts, fraction = TRUE, aggregated = FALSE, notDetected = "", sep = ";") {

  # get metadata and ident class
  if (class(counts)[1] == "Seurat") {
    counts <- counts@meta.data
    aggregated <- TRUE
  } else if (class(counts)[1] == "SingleCellExperiment") {
    counts <- as.data.frame(counts@colData)
    aggregated <- TRUE
  }

  if (aggregated) {
    # in case barcodes were previously already aggregated per cell
    lineagePerCell.dist.df <- counts %>%
      tibble::rownames_to_column("cellid") %>%
      dplyr::select(.data$cellid, .data$barcode) %>%
      dplyr::mutate(barcode = ifelse(.data$barcode == notDetected, NA, .data$barcode)) %>%
      dplyr::mutate(number_of_lineage_barcodes = stringr::str_count(.data$barcode, pattern = sep) + 1) %>%
      dplyr::mutate(number_of_lineage_barcodes = ifelse(is.na(.data$number_of_lineage_barcodes), 0, .data$number_of_lineage_barcodes)) %>%
      dplyr::select(-.data$barcode)
  } else {
    # Detected lineage barcodes per cell
    lineagePerCell.dist.df <- counts %>%
      dplyr::select(.data$cellid, .data$barcode) %>%
      dplyr::group_by(.data$cellid) %>%
      dplyr::tally(., name = "number_of_lineage_barcodes")
  }

  lineagePerCell.dist.df <- lineagePerCell.dist.df %>%
    dplyr::count(.data$number_of_lineage_barcodes) %>%
    dplyr::mutate(frac = .data$n / nrow(lineagePerCell.dist.df))

  if (fraction) {
    p <- lineagePerCell.dist.df %>%
      ggplot2::ggplot(aes(x = .data$number_of_lineage_barcodes, y = .data$frac))
  } else {
    p <- lineagePerCell.dist.df %>%
      ggplot2::ggplot(aes(x = .data$number_of_lineage_barcodes, y =.data$n))
  }
  p <- p +
    ggplot2::theme(
      axis.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold")
    ) +
    ggplot2::geom_bar(stat = "identity",
                      fill = "blue",
                      width = .75) +
    ggplot2::xlab("# Barcodes") +
    ggplot2::ggtitle(paste("Number of barcodes detected per cell")) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = integer_breaks())

  if (fraction) {
    p <- p +
      ylab("Fraction of cells")
  } else {
    p <- p +
      ylab("# cells")
  }

  return(p)
}


#' plotUmiPerBarcode
#'
#' Plot number of UMIs supporting the most frequenc barcode per cell
#'
#' @param counts Dataframe with barcodes and UMI counts per cell
#' @param fraction Boolean, whether to plot fraction or number of cells
#'
#' @return Returns a plot
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export

plotUmiPerBarcode <- function(counts, fraction = TRUE) {
  # Number of UMIs supporting the most frequent barcode
  max.umi.per.cell <- counts %>%
    dplyr::group_by(.data$cellid) %>%
    dplyr::summarise(max = max(.data$bc.umi.count))

  max.umi.per.cell <- max.umi.per.cell %>%
    dplyr::count(.data$max) %>%
    dplyr::mutate(frac = .data$n / nrow(max.umi.per.cell))

  if (fraction) {
    p <- max.umi.per.cell %>%
      ggplot2::ggplot(aes(x = .data$max, y = .data$frac))
  } else {
    p <- max.umi.per.cell %>%
      ggplot2::ggplot(aes(x = .data$max, y = .data$n))
  }
  p <- p +
    ggplot2::theme(
      axis.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold")
    ) +
    ggplot2::geom_bar(stat = "identity", width = .75) +
    ggplot2::xlab("# UMI") +
    ggplot2::ggtitle("Number of UMI supporting the most frequent barcode") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = integer_breaks())

  if (fraction) {
    p <- p +
      ylab("Fraction of cells")
  } else {
    p <- p +
      ylab("# cells")
  }

  return(p)
}


#' plotUmiFilterThresholds
#'
#' Plot number of cells with barcodes annotated given different UMI count filter thresholds
#'
#' @param counts Dataframe with barcodes and UMI counts per cell
#'
#' @return Returns a plot
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export
#'
plotUmiFilterThresholds <- function(counts) {
  max.umi.per.cell <- counts %>%
    # get UMI count of dominant barcode per cell
    dplyr::group_by(.data$cellid) %>%
    dplyr::summarise(max = max(.data$bc.umi.count)) %>%
    # cummulative count of cells
    dplyr::group_by(.data$max) %>%
    dplyr::summarise(count_max = dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(.data$max)) %>%
    dplyr::mutate(cumsum_count_max = cumsum(.data$count_max))

  p <-
    ggplot2::ggplot(max.umi.per.cell, aes(x = .data$max, y = .data$cumsum_count_max)) +
    ggplot2::geom_bar(stat = "identity", width = .75) +
    ggplot2::ylab("# cells") +
    ggplot2::xlab("UMI threshold") +
    ggplot2::ggtitle("Number of cells with barcodes annotated given UMI threshold") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = integer_breaks())

  return(p)
}
