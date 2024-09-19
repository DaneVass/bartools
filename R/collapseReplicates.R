#' @title
#' Collapse technical replicates
#'
#' @description
#' Collapse technical replicates in a DGEList object by mean or sum. Updates lib.size to the read sum per collapsed sample.
#' Modified from collapseReplicates function from DESeq2 package to accept DGEList objects and to allow collapsing replicated by mean or sum.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param group Column name in sample metadata to group samples by (string).
#' @param sampleNames Optional, column name in sample metadata that contains unique sample names (string). If provided, a new column collapsed_samples will be added to the sample metadata with the collapsed samples
#' @param renameCols Whether to rename the columns of the returned object using the levels of the grouping factor (boolean). Default = `TRUE`.
#' @param showReps Whether to print replicate column names to console (boolean). Default = `FALSE`.
#' @param method Method to collapse replicates by, one of `mean` or `sum` (string). Default = `mean`.
#'
#' @return Returns a DGE object with normalised counts for a sample.
#' @export
#'
#' @examples
#' data(test.dge)
#' collapseReplicates(test.dge, group = "group", method = "mean")
#'

collapseReplicates <-
  function (dgeObject,
            group,
            sampleNames = NULL,
            renameCols = TRUE,
            showReps = FALSE,
            method = "mean") {
    inputChecks(dgeObject, groups = c(group, sampleNames))

    counts <- dgeObject$counts

    # check inputs
    if (!method %in% c("mean", "sum")) {
      stop("Choose either mean or sum as method to collapse replicates.")
    }

    group <- dgeObject$samples[[group]]
    if (!is.null(sampleNames)) {
      sampleNames <- dgeObject$samples[[sampleNames]]
    }

    if (!is.factor(group)) {
      group <- factor(group)
    }
    group <- droplevels(group)

    sp <- split(seq(along = group), group)
    if (isTRUE(showReps)) {
      print(sp)
    }

    # combine by mean or sum
    if (method == "mean") {
      counts <-
        sapply(sp, function(i)
          rowMeans(counts[, i, drop = FALSE]))
    }
    if (method == "sum") {
      counts <-
        sapply(sp, function(i)
          rowSums(counts[, i, drop = FALSE]))
    }

    # select columns to keep
    # this selects always the first column in each group.
    # this assumes that all values are the same across samples in a group.
    colsToKeep <- sapply(sp, `[`, 1)
    collapsed <- dgeObject[, colsToKeep]
    dimnames(counts) <- dimnames(collapsed)
    collapsed$counts <- counts

    if (!is.null(sampleNames)) {
      collapsed$samples$collapsed_samples <-
        sapply(sp, function(i)
          paste(sampleNames[i],
                collapse = ","))
    }

    # rename cols true by default
    if (renameCols) {
      colnames(collapsed) <- levels(group)
    }

    # update lib.size
    # important when using normalization functions that use this column, like edgeR::cpm()
    collapsed$samples$lib.size <- colSums(collapsed$counts)

    return(collapsed)
  }
