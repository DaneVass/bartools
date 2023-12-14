#' Calculate correlation between technical replicates
#'
#' Calculate and return list of correlation between paired technical replicates in a dataset
#'
#' @param dgeObject DGEList object with barcode counts with technical replicates.
#' @param group Column name in sample metadata with replicate information (string). Correlation can only be calculated between pairs of replicates.
#' @param corrThreshold Threshold distinguishing good vs bad correlation between technical replicates (decimal). Default = `0.8`.
#' @param return Which values to return, one of `good`, `bad`, `all` (string). Default = `all`.
#' @param method Method for correlation, one of `spearman`, `pearson` (string). Default = `pearson`.
#' @param ignoreZero Remove barcodes where both replicates have 0 counts before calculating correlation (boolean). Including zeros can increase the spearman correlation, especially for samples with few barcodes detected. Default = `FALSE`.
#'
#' @return Returns a plot of the read counts per barcode (row) in a data frame
#' @export
#'
#' @examples
#' data(test.dge)
#' calcReplicateCorr(test.dge, group = "group")

# get paired technical replicates
calcReplicateCorr <-
  function(dgeObject,
           group = NULL,
           corrThreshold = 0.8,
           return = "all",
           method = "pearson",
           ignoreZero = FALSE) {
    inputChecks(dgeObject, groups = group)

    if (!return %in% c("all", "good", "bad")) {
      stop("Return must be one of all, good or bad")
    }

    counts <- dgeObject$counts

    if (is.null(group)) {
      stop(
        "Please supply column name in dgeObject$samples dataframe with replicate information"
      )
    }

    # define grouping column in dge sample metadata
    singletons <- which(table(dgeObject$samples[, group]) == 1)
    if (length(singletons) > 0) {
      message(paste0(
        "Following groups have only one technical replicate: ",
        paste(names(singletons), collapse = ", ")
      ))
    }
    multiple_replicates <-
      which(table(dgeObject$samples[, group]) > 2)
    if (length(multiple_replicates) > 0) {
      message(paste0(
        "Following groups have more than 2 replicates: ",
        paste(names(multiple_replicates), collapse = ", ")
      ))
    }
    paired <- which(table(dgeObject$samples[, group]) == 2)

    # filter samples with technical replicates
    keep <- which(dgeObject$samples[, group] %in% (names(paired)))
    dgeObject.filter <- dgeObject[, keep]

    corrs <- lapply(unique(names(paired)), function(x) {
      data <-
        as.data.frame(dgeObject$counts[, dgeObject$samples[, group] == as.character(x)])

      if (ignoreZero == TRUE) {
        data <- data[rowSums(data) != 0,]
      }
      # # Previous implementation, more in line with plotBarcodeRegression
      # # But not clear why adjusted r is needed. Also returns strong negative correlation as positive value.
      # fit <- stats::lm(data[, 1] ~ data[, 2])
      # adj.r2 <- summary(fit)$adj.r.squared
      # corr <- sqrt(adj.r2)

      return(cor(data[, 1], data[, 2], method = method))
    })

    names(corrs) <- unique(names(paired))
    corrs <- unlist(corrs)

    # good corrs
    good.corrs <- corrs[which(corrs > corrThreshold)]
    # bad corrs
    bad.corrs <- corrs[which(corrs < corrThreshold)]

    if (return == "all") {
      return(corrs)
    }
    if (return == "bad") {
      return(bad.corrs)
    }
    if (return == "good") {
      return(good.corrs)
    }
  }
