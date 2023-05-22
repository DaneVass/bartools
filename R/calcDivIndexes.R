#' calcDivIndexes
#'
#' Calculates common diversity indices per sample in a dataframe of counts
#'
#' @title
#' Calculate diversity indices for barcode samples
#'
#' @description
#' Takes a dataframe of barcode counts and computes shannon,
#' simpson, inverse simpson and gini coefficients for each
#' sample
#'
#' @param counts dataframe of raw or normalised counts with samples as columns and observations/barcodes as rows
#' @param group vector of sample group metadata, assumes the same sample order as counts
#'
#' @return Returns a data-frame containing the calculated diversity indices per sample
#' @importFrom vegan diversity
#' @importFrom ineq Gini
#' @export
#' @examples
#' indexes <- calcDivIndexes(counts = test.dge$counts)

calcDivIndexes <- function(counts, group = NULL){
  colnames <- colnames(counts)
  counts <- as.data.frame(counts)
  df <- data.frame(name = c(), shannon = as.numeric(), simpson = as.numeric(),
                     invsimpson = as.numeric(), gini = as.numeric())

  for (i in 1:ncol(counts)) {
    name <- colnames(counts)[i]
    shannon <- vegan::diversity(counts[, i], index = "shannon")
    simpson <- vegan::diversity(counts[, i], index = "simpson")
    invsimpson <- vegan::diversity(counts[, i], index = "invsimpson")
    gini <- ineq::Gini(counts[, i])
    x <- data.frame(name, shannon, simpson, invsimpson, gini)
    df <- rbind(df, x)
  }
  #print(df)
  
  # add group metadata
  if (is.null(group)){
    group <- rep("Group1", ncol(counts))
    Group <- as.data.frame(colnames(counts))
    Group$group <- group
    colnames(Group)[1] <- "name"
  } else {
    Group <- as.data.frame(colnames(counts))
    Group$group <- group
    colnames(Group)[1] <- "name"
  }
  df <- dplyr::left_join(df, Group)
  df$name <- as.factor(df$name)
  df$group <- as.factor(df$group)
  return(df)
}
