#' plotBARtabFilterQC
#'
#' Plot filtered read percentages from a BARtab run folder
#'
#' @param dir directory where BARtab was successfuly run on barcode count datasets
#' @param pattern.log regex string to specify filter stage log files
#' @param pattern.value grep on this string in log files
#' @param full.names Logical. Return full names of files detected by regex search
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
                               pattern.log = "*filter.log",
                               pattern.value = "reads",
                               full.names = T,
                               normalised = F,
                               plot = T,
                               title = "BARtab Filter QC"){

  # get log files from directory
  logs <- list.files(dir, pattern = pattern.log, full.names = full.names, recursive = recursive)

  # setup output data.frame to store filter information
  final.df <- data.frame(sample = as.character(), input = as.numeric(), output = as.numeric(), discarded = as.numeric(), pct.kept = as.numeric())

  # loop over log files and pull out relevant filtering stats
  for (log in logs){

    # get samplename
    samp <- utils::tail(stringr::str_split(log, pattern = "/")[[1]], 1)
    samp <- utils::head(stringr::str_split(samp, pattern = "\\.")[[1]][1],1)

    # get values
    filtered.reads <- grep(pattern.value, readLines(log), value = T)
    filtered.uniq <- stringr::str_extract(filtered.reads, '[0-9]{3,20}')
    input <- as.numeric(filtered.uniq[1])
    output <- as.numeric(filtered.uniq[2])
    discarded <- as.numeric(filtered.uniq[3])
    pct.kept <- 100*(output/input)
    # concatenate to dataframe
    df <- data.frame(sample = samp, input = input, output = output, discarded = discarded, pct.kept = pct.kept)
    final.df <- rbind(final.df, df)
  }

  # factorise samplenames
  final.df$sample <- as.factor(gsub("-", "_",final.df$sample))

  if (plot){
    # plot data if requested
    plot.dat <- final.df
    plot.dat$pct.kept <- NULL
    plot.dat$output <- NULL
    plot.dat <- reshape2::melt(plot.dat)

    # normalise if requested
    if(normalised){
      p <- ggplot(plot.dat, aes(fill=variable, y=value, x=sample)) +
        geom_bar(position="fill", stat="identity") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = c("dodgerblue2", "lightblue")) +
        ggplot2::ggtitle("Percentage of barcode reads kept post filtering")
      print(p)
    } else {
      p <- ggplot(plot.dat, aes(fill=variable, y=value, x=sample)) +
        geom_bar(position="stack", stat="identity") +
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
#' @param pattern.log regex string to specify filter stage log files
#' @param full.names Logical. Return full names of files detected by regex search
#' @param recursive Logical. TRUE will recurse regex search into subdirectories#'
#' @param plot Logical. Generate plots. False returns raw data
#' @param title Optional. title of plots.
#'
#' @return Returns a plot / dataset containing BARtab Filtering stage QC data
#'
#' @export

plotBARtabMapQC <- function(dir = NULL,
                            recursive = T,
                            pattern.log = "*bowtie.log",
                            full.names = T,
                            plot = T,
                            title = "BARtab Mapping QC"){

  # get log files from directory
  logs <- list.files(dir, pattern = pattern.log, full.names = full.names, recursive = recursive)

  # setup output data.frame to store filter information
  final.df <- data.frame(sample = as.character(), percent = as.numeric())

  for (log in logs){

    # get samplename from log
    samp <- tail(stringr::str_split(log, pattern = "/")[[1]],1)
    samp <- stringr::str_split(samp, pattern = "_", n = 2)[[1]][1]

    # filter alignment info from logfile
    aligned.list <- list(grep("#", readLines(log), value = T))
    aligned.uniq <- aligned.list[[1]][2]
    aligned.uniq <- stringr::str_extract(aligned.uniq, '\\([0-9][0-9].[0-9][0-9]%\\)')
    aligned.uniq <- gsub("\\(", "", aligned.uniq)
    aligned.uniq <- gsub("\\)", "", aligned.uniq)
    aligned.uniq <- gsub("%", "", aligned.uniq)
    df <- data.frame(sample = samp, percent = as.numeric(aligned.uniq))
    final.df <- rbind(final.df, df)
  }

  # factorise sample group
  final.df$sample <- as.factor(gsub("-", "_",final.df$sample))
  final.df$group <- as.factor(gsub("_PCR[12]$", "",final.df$sample))
  final.df$group <- as.factor(gsub("_[12]$", "",final.df$sample))

  colour_breaks <- c(0, 20, 40, 60, 80, 100)
  colours <- rev(c("darkblue", "lightblue", "yellow", "orange", "red", "firebrick"))
  
  # plot data
  if (plot) {
    p <- ggplot2::ggplot(final.df, ggplot2::aes(sample, percent, fill = percent)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Percentage of barcode reads aligning to the reference library") +
      #ggplot2::theme(legend.position = "none") +
      scale_fill_gradientn(
        limits  = range(0,100),
        colours = colours[c(1, seq_along(colours), length(colours))],
        values  = c(0, scales::rescale(colour_breaks, from = range(0,100)), 1),
      )
    
    suppressWarnings(print(p))
  } else {
    return(final.df)
  }
}

