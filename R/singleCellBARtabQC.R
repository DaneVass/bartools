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
  counts <- read.delim(countsPath, header = TRUE)
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
#' Aggregate barcodes and barcode UMIs per cell
#'
#' @param counts Dataframe with barcodes and UMI counts per cell
#' @param sep Separator between barcodes (string). Default = `;`.
#'
#' @return Returns a data frame with one row per cell ID
#'
#' @export
#
aggregateBarcodes <- function(counts, sep = ";") {
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
#'
#' @export
#
filterBarcodes <-
  function(counts,
           umiCountFilter = 1,
           umiFractionFilter = 0.3) {
    counts_filtered <- counts %>%
      dplyr::filter(bc.umi.count >= umiCountFilter)

    counts_filtered <- counts_filtered %>%
      dplyr::group_by(cellid) %>%
      dplyr::mutate(max_count = max(bc.umi.count)) %>%
      dplyr::arrange(cellid) %>%
      dplyr::filter(bc.umi.count / umiFractionFilter >= max_count)

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
#'
#' @export

plotBarcodesPerCell <- function(counts, fraction = TRUE, aggregated = FALSE, notDetected = "", sep = ";") {

  # get metadata and ident class
  if (class(counts)[1] == "Seurat") {
    counts <- counts@meta.data
    aggregated <- TRUE
  } else if (class(counts)[1] == "SingleCellExperiment") {
    counts <- counts@colData
    aggregated <- TRUE
  }

  if (aggregated) {
    # in case barcodes were previously already aggregated per cell
    lineagePerCell.dist.df <- counts %>%
      tibble::rownames_to_column("cellid") %>%
      dplyr::select(cellid, barcode) %>%
      dplyr::mutate(barcode = ifelse(barcode == notDetected, NA, barcode)) %>%
      dplyr::mutate(number_of_lineage_barcodes = stringr::str_count(barcode, pattern = sep) + 1) %>%
      dplyr::mutate(number_of_lineage_barcodes = ifelse(is.na(number_of_lineage_barcodes), 0, number_of_lineage_barcodes)) %>%
      dplyr::select(-barcode)
  } else {
    # Detected lineage barcodes per cell
    lineagePerCell.dist.df <- counts %>%
      dplyr::select(cellid, barcode) %>%
      dplyr::group_by(cellid) %>%
      dplyr::tally(., name = "number_of_lineage_barcodes")
  }

  lineagePerCell.dist.df <- lineagePerCell.dist.df %>%
    dplyr::count(number_of_lineage_barcodes) %>%
    dplyr::mutate(frac = n / nrow(lineagePerCell.dist.df))

  if (fraction) {
    p <- lineagePerCell.dist.df %>%
      ggplot2::ggplot(aes(x = number_of_lineage_barcodes, y = frac))
  } else {
    p <- lineagePerCell.dist.df %>%
      ggplot2::ggplot(aes(x = number_of_lineage_barcodes, y = n))
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
#'
#' @export

plotUmiPerBarcode <- function(counts, fraction = TRUE) {
  # Number of UMIs supporting the most frequent barcode
  max.umi.per.cell <- counts %>%
    dplyr::group_by(cellid) %>%
    dplyr::summarise(max = max(bc.umi.count))

  max.umi.per.cell <- max.umi.per.cell %>%
    dplyr::count(max) %>%
    dplyr::mutate(frac = n / nrow(max.umi.per.cell))

  if (fraction) {
    p <- max.umi.per.cell %>%
      ggplot2::ggplot(aes(x = max, y = frac))
  } else {
    p <- max.umi.per.cell %>%
      ggplot2::ggplot(aes(x = max, y = n))
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
#'
#' @export
#'
plotUmiFilterThresholds <- function(counts) {
  max.umi.per.cell <- counts %>%
    # get UMI count of dominant barcode per cell
    dplyr::group_by(cellid) %>%
    dplyr::summarise(max = max(bc.umi.count)) %>%
    # cummulative count of cells
    dplyr::group_by(max) %>%
    dplyr::summarise(count_max = dplyr::n()) %>%
    dplyr::arrange(desc(max)) %>%
    dplyr::mutate(cumsum_count_max = cumsum(count_max))

  p <-
    ggplot2::ggplot(max.umi.per.cell, aes(x = max, y = cumsum_count_max)) +
    ggplot2::geom_bar(stat = "identity", width = .75) +
    ggplot2::ylab("# cells") +
    ggplot2::xlab("UMI threshold") +
    ggplot2::ggtitle("Number of cells with barcodes annotated given UMI threshold") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = integer_breaks())

  return(p)
}
