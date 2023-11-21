#' calcDivIndexes
#'
#' Calculates common diversity indices per sample in a dataframe of counts
#'
#' @title
#' Calculate diversity indices for barcode samples
#'
#' @description
#' Takes a DGEList object of barcode counts and computes shannon,
#' simpson, inverse simpson and gini coefficients for each
#' sample
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param group Optional, column name in sample metadata to group samples by (string).
#'
#' @return Returns a data-frame containing the calculated diversity indices per sample
#' @importFrom vegan diversity
#' @importFrom ineq Gini
#' @export
#' @examples
#' data(test.dge)
#' indexes <- calcDivIndexes(test.dge)

calcDivIndexes <- function(dgeObject, group = NULL) {
  inputChecks(dgeObject, groups = group)

  counts <- as.data.frame(dgeObject$counts)

  df <-
    data.frame(
      name = c(),
      shannon = as.numeric(),
      simpson = as.numeric(),
      invsimpson = as.numeric(),
      gini = as.numeric()
    )

  for (i in 1:ncol(counts)) {
    name <- colnames(counts)[i]
    shannon <- vegan::diversity(counts[, i], index = "shannon")
    simpson <- vegan::diversity(counts[, i], index = "simpson")
    invsimpson <-
      vegan::diversity(counts[, i], index = "invsimpson")
    gini <- ineq::Gini(counts[, i])
    x <- data.frame(name, shannon, simpson, invsimpson, gini)
    df <- rbind(df, x)
  }

  # add group metadata
  if (is.null(group)) {
    group <- rep("Group1", ncol(counts))
    Group <- as.data.frame(colnames(counts))
    Group$group <- group
    colnames(Group)[1] <- "name"
  } else {
    Group <- as.data.frame(colnames(counts))
    Group$group <- dgeObject$samples[group]
    colnames(Group)[1] <- "name"
  }
  df <- dplyr::left_join(df, Group, by = "name")
  df$name <- as.factor(df$name)
  df$group <- as.factor(df$group)
  return(df)
}
