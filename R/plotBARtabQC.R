

#' plotBARtabFilterQC
#'
#' Plot filtered read percentages from a BARtab run folder
#'
#' @param dir directory where BARtab was successfuly run on barcode count datasets
#' @param patternLog regex string to specify filter stage log files
#' @param fullNames Logical. Return full names of files detected by regex search
#' @param recursive Logical. TRUE will recurse regex search into subdirectories
#' @param normalised Logical. log10 normalise counts
#' @param plot Logical. Generate plots. False returns raw data
#' @param title Optional. title of plots.
#'
#' @return Returns a plot / dataset containing BARtab Filtering stage QC data
#'
#' @importFrom rlang .data
#' @export

plotBARtabFilterQC <- function(dir = NULL,
                               recursive = T,
                               patternLog = "*filter.log",
                               fullNames = T,
                               normalised = F,
                               plot = T,
                               title = "BARtab Filter QC") {
  # get log files from directory
  logs <-
    list.files(dir,
               pattern = patternLog,
               full.names = fullNames,
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
    log_lines <- readLines(log)
    # make it work with old fastx-toolkit log and new fastp log
    input_reads <- grep("(Input: [0-9]* reads.|total reads: [0-9*])", log_lines, value = T)[1]
    input_reads <- as.numeric(stringr::str_extract(input_reads, pattern = '\\d+'))

    passed_reads <- grep("(Output: [0-9]* reads.|reads passed filter: [0-9*])", log_lines, value = T)[1]
    passed_reads <- as.numeric(stringr::str_extract(passed_reads, pattern = '\\d+'))

    discarded_reads <- input_reads - passed_reads
    pct.kept <- 100 * (passed_reads / input_reads)
    # concatenate to dataframe
    df <-
      data.frame(
        sample = samp,
        input = input_reads,
        output = passed_reads,
        discarded = discarded_reads,
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
      p <- ggplot2::ggplot(plot.dat, ggplot2::aes(fill = .data$variable, y = .data$value, x = .data$sample)) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = c("dodgerblue2", "lightblue")) +
        ggplot2::ggtitle("Percentage of barcode reads kept post filtering")
      return(p)
    } else {
      p <- ggplot2::ggplot(plot.dat, ggplot2::aes(fill = .data$variable, y = .data$value, x = .data$sample)) +
        ggplot2::geom_bar(position = "stack", stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = c("dodgerblue2", "lightblue")) +
        ggplot2::ggtitle("Total barcode reads kept post filtering")
      return(p)
    }
  } else {
    return(final.df)
  }
}

#' plotBARtabMapQC
#'
#' Plot mapped read percentages from a BARtab run folder. Shows reads with at least one alignment.
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
#' @importFrom rlang .data
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
               full.names = fullNames,
               recursive = recursive)

  # setup output data.frame to store filter information
  final.df <-
    data.frame(sample = as.character(), percent = as.numeric())

  for (log in logs) {
    # get samplename from log
    samp <- utils::tail(stringr::str_split(log, pattern = "/")[[1]], 1)
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
      ggplot2::ggplot(final.df, ggplot2::aes(.data$sample, .data$percent, fill = .data$percent)) +
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
