

# function to check if general inputs are valid
inputChecks <-
  function(dgeObject,
           barcodes = NULL,
           samples = NULL,
           groups = NULL,
           conditions = NULL) {
    # check that input is DGEList object
    if (methods::is(dgeObject)[1] != "DGEList") {
      stop("Please supply a valid DGEList object as input")
    }

    # check if selected barcodes are present
    if (!is.null(barcodes)) {
      if (!all(barcodes %in% rownames(dgeObject))) {
        stop(paste(
          "Following barcodes are not present in the input DGEList object:",
          paste0(barcodes[!barcodes %in% rownames(dgeObject)], collapse = ", ")
        ))
      }
    }

    # check if selected samples are present
    if (!is.null(samples)) {
      if (!all(samples %in% colnames(dgeObject))) {
        stop(paste(
          "Following samples are not in input data:",
          paste0(samples[!samples %in% colnames(dgeObject)], collapse = ", ")
        ))
      }
    }

    # check if group is in samples metadata
    if (!is.null(groups)) {
      if (!all(groups %in% colnames(dgeObject$samples))) {
        stop(paste(
          "Following groups are not a column in the samples metadata:",
          paste0(groups[!groups %in% colnames(dgeObject$samples)], collapse = ", ")
        ))
      }
    }

    # check if selected conditions are in group column
    if (!is.null(conditions)) {
      if (is.null(groups)) {
        stop("Specify column to group by if selecting conditions.")
      }
      if (!all(conditions %in% dgeObject$samples[[groups]])) {
        stop(paste(
          "Following conditions are not in group column:",
          paste0(conditions[!conditions %in% dgeObject$samples[[groups]]], collapse = ", ")
        ))
      }
    }
  }

#### 2 functions used by plotBarcodeCorrelation and plotBarcodeDistance
#' Cluster correlation matrix
#'
#' cluster a correlation matrix using hierarchical clustering
#'
#' @param cormat matrix of correlation values
#'
#' @return Returns a matrix of correlation values with columns and rows hierarchically clustered
#'
#' @export
#'
cluster_cormat <- function(correlation_mat) {
  dd <- stats::as.dist((1 - correlation_mat) / 2)
  hc <- stats::hclust(dd)
  correlation_mat <- correlation_mat[hc$order, hc$order]
  return(correlation_mat)
}

#' Get upper triangle
#'
#' change lower triangle values in a mirrored matrix to NA
#'
#' @param correlation_mat matrix of correlation values
#'
#' @return Returns a matrix of correlation values with lower triangle values changed to NA
#' @export
#'
get_upper_tri <- function(correlation_mat) {
  correlation_mat[lower.tri(correlation_mat)] <- NA
  return(correlation_mat)
}
