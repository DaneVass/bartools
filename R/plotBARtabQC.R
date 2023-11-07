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
#'
#' @return Returns a data frame with one row per cell ID
#'
#' @export
#
aggregateBarcodes <- function(counts) {
  bc.counts <- as.data.frame(data.table::dcast(
    data.table::setDT(counts),
    cellid ~ .,
    value.var = c("barcode", "bc.umi.count"),
    fun.aggregate = function(x) {
      paste(x, collapse = ";")
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
#' Plot number of detected barcodes per cell
#'
#' @param counts Dataframe with barcodes and UMI counts per cell
#' @param fraction Boolean, whether to plot fraction or number of cells
#'
#' @return Returns a plot
#'
#' @export

plotBarcodesPerCell <- function(counts, fraction = TRUE) {
  # Detected lineage barcodes per cell
  lineagePerCell.dist.df <- counts %>%
    dplyr::select(cellid, barcode) %>%
    dplyr::group_by(cellid) %>%
    dplyr::tally(., name = "number_of_lineage_barcodes") %>%
    dplyr::arrange(dplyr::desc(number_of_lineage_barcodes))

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


#' plotBARtabFilterQC
#'
#' Plot filtered read percentages from a BARtab run folder
#'
#' @param dir directory where BARtab was successfuly run on barcode count datasets
#' @param patternLog regex string to specify filter stage log files
#' @param pattern.value grep on this string in log files
#' @param fullNames Logical. Return full names of files detected by regex search
#' @param recursive Logical. TRUE will recurse regex search into subdirectories#'
#' @param normalised Logical. log10 normalise counts
#' @param plot Logical. Generate plots. False returns raw data
#' @param title Optional. title of plots.
#'
#' @return Returns a plot / dataset containing BARtab Filtering stage QC data
#'
#' @export

plotBARtabFilterQC <- function(dir = NULL,
                               recursive = T,
                               patternLog = "*filter.log",
                               pattern.value = "reads",
                               fullNames = T,
                               normalised = F,
                               plot = T,
                               title = "BARtab Filter QC") {
  # get log files from directory
  logs <-
    list.files(dir,
               pattern = patternLog,
               fullNames = fullNames,
               recursive = recursive)

  # setup output data.frame to store filter information
  final.df <-
    data.frame(
      sample = as.character(),
      input = as.numeric(),
      output = as.numeric(),
      discarded = as.numeric(),
      pct.kept = as.numeric()
    )

  # loop over log files and pull out relevant filtering stats
  for (log in logs) {
    # get samplename
    samp <-
      utils::tail(stringr::str_split(log, pattern = "/")[[1]], 1)
    samp <-
      utils::head(stringr::str_split(samp, pattern = "\\.")[[1]][1], 1)

    # get values
    filtered.reads <- grep(pattern.value, readLines(log), value = T)
    filtered.uniq <-
      stringr::str_extract(filtered.reads, '[0-9]{3,20}')
    input <- as.numeric(filtered.uniq[1])
    output <- as.numeric(filtered.uniq[2])
    discarded <- as.numeric(filtered.uniq[3])
    pct.kept <- 100 * (output / input)
    # concatenate to dataframe
    df <-
      data.frame(
        sample = samp,
        input = input,
        output = output,
        discarded = discarded,
        pct.kept = pct.kept
      )
    final.df <- rbind(final.df, df)
  }

  # factorise samplenames
  final.df$sample <- as.factor(gsub("-", "_", final.df$sample))

  if (plot) {
    # plot data if requested
    plot.dat <- final.df
    plot.dat$pct.kept <- NULL
    plot.dat$output <- NULL
    plot.dat <- reshape2::melt(plot.dat)

    # normalise if requested
    if (normalised) {
      p <- ggplot(plot.dat, aes(fill = variable, y = value, x = sample)) +
        geom_bar(position = "fill", stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = c("dodgerblue2", "lightblue")) +
        ggplot2::ggtitle("Percentage of barcode reads kept post filtering")
      print(p)
    } else {
      p <- ggplot(plot.dat, aes(fill = variable, y = value, x = sample)) +
        geom_bar(position = "stack", stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = c("dodgerblue2", "lightblue")) +
        ggplot2::ggtitle("Total barcode reads kept post filtering")
      print(p)
    }
  } else {
    return(final.df)
  }
}

#' plotBARtabMapQC
#'
#' Plot mapped read percentages from a BARtab run folder
#'
#' @param dir directory where BARtab was successfuly run on barcode count datasets
#' @param patternLog regex string to specify filter stage log files
#' @param fullNames Logical. Return full names of files detected by regex search
#' @param recursive Logical. TRUE will recurse regex search into subdirectories#'
#' @param plot Logical. Generate plots. False returns raw data
#' @param title Optional. title of plots.
#'
#' @return Returns a plot / dataset containing BARtab Filtering stage QC data
#'
#' @export

plotBARtabMapQC <- function(dir = NULL,
                            recursive = T,
                            patternLog = "*bowtie.log",
                            fullNames = T,
                            plot = T,
                            title = "BARtab Mapping QC") {
  # get log files from directory
  logs <-
    list.files(dir,
               pattern = patternLog,
               fullNames = fullNames,
               recursive = recursive)

  # setup output data.frame to store filter information
  final.df <-
    data.frame(sample = as.character(), percent = as.numeric())

  for (log in logs) {
    # get samplename from log
    samp <- tail(stringr::str_split(log, pattern = "/")[[1]], 1)
    samp <- stringr::str_split(samp, pattern = "_", n = 2)[[1]][1]

    # filter alignment info from logfile
    aligned.list <- list(grep("#", readLines(log), value = T))
    aligned.uniq <- aligned.list[[1]][2]
    aligned.uniq <-
      stringr::str_extract(aligned.uniq, '\\([0-9][0-9].[0-9][0-9]%\\)')
    aligned.uniq <- gsub("\\(", "", aligned.uniq)
    aligned.uniq <- gsub("\\)", "", aligned.uniq)
    aligned.uniq <- gsub("%", "", aligned.uniq)
    df <-
      data.frame(sample = samp, percent = as.numeric(aligned.uniq))
    final.df <- rbind(final.df, df)
  }

  # factorise sample group
  final.df$sample <- as.factor(gsub("-", "_", final.df$sample))
  final.df$group <- as.factor(gsub("_PCR[12]$", "", final.df$sample))
  final.df$group <- as.factor(gsub("_[12]$", "", final.df$sample))

  colour_breaks <- c(0, 20, 40, 60, 80, 100)
  colours <- c("darkblue", "lightblue", "yellow", "orange", "red2")

  # plot data
  if (plot) {
    p <-
      ggplot2::ggplot(final.df, ggplot2::aes(sample, percent, fill = percent)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Percentage of barcode reads aligning to the reference library") +
      #ggplot2::theme(legend.position = "none") +
      scale_fill_gradientn(
        limits  = range(0, 100),
        colours = colours[c(1, seq_along(colours), length(colours))],
        values  = c(0, scales::rescale(colour_breaks, from = range(0, 100)), 1),
      )

    suppressWarnings(print(p))
  } else {
    return(final.df)
  }
}
