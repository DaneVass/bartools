#' plotSampleCumSum
#'
#' @title
#' Plot the cumulative sum of barcode abundances within samples to visualize sample diversity.
#'
#' @param dgeObject DGE object with barcode counts
#' @param samples Sample names
#'
#' @return Returns a cumulative sum plot
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export

plotSampleCumSum <- function(dgeObject,
                             samples = NULL) {
  inputChecks(dgeObject, samples = samples)

  counts <- dgeObject$counts
  if (!is.null(samples)) {
    counts <- counts[, samples]
  }

  counts <- counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Barcode") %>%
    tidyr::pivot_longer(-.data$Barcode, names_to = "sample", values_to = "count")

  counts <- counts %>%
    dplyr::arrange(dplyr::desc(.data$count)) %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::mutate(freq = .data$count / sum(.data$count)) %>%
    # remove all samples not detected within sample
    dplyr::filter(.data$freq > 0) %>%
    dplyr::mutate(rank = dplyr::row_number(), cumsum = cumsum(.data$freq))

  # shuffle points so samples are not hidden
  counts[sample(1:nrow(counts)),]  %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$rank, y = .data$cumsum, color = .data$sample)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::ylab("Cumulative Sum (Proportion)") +
    ggplot2::xlab("Barcode rank within sample") +
    ggplot2::theme_bw()
}
