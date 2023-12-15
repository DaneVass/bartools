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
    tidyr::pivot_longer(-Barcode, names_to = "sample", values_to = "count")

  counts <- counts %>%
    arrange(dplyr::desc(count)) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(freq = count / sum(count)) %>%
    # remove all samples not detected within sample
    dplyr::filter(freq > 0) %>%
    dplyr::mutate(rank = dplyr::row_number(), cumsum = cumsum(freq))

  # shuffle points so samples are not hidden
  counts[sample(1:nrow(counts)),]  %>%
    ggplot(aes(x = rank, y = cumsum, color = sample)) +
    geom_point(alpha = 0.5) +
    ylab("Cumulative Sum (Proportion)") +
    xlab("Barcode rank within sample") +
    theme_bw()
}
