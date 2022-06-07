#' Calculate correlation between technical replicates
#'
#' Calculate and return list of correlation between paired technical replicates in a dataset
#'
#' @param dge DGEList object containing technical replicates to be compared
#' @param group grouping variable from metadata containing technical replicate information
#' @param corr.thresh threshold distinguishing good vs bad correlation between technical replicates
#' @param return which values to return. One of "good", "bad", "all"
#'
#' @return Returns a plot of the read counts per barcode (row) in a data frame
#' @export
#'
#' @examples
#' data(test.dge)
#' calcReplicateCorr(test.dge, group = "group")

# get paired technical replicates
calcReplicateCorr <- function(dge, group, corr.thresh = 0.8, return = "all"){
  #if(return != "all" | return != "good" | return != "bad"){
  #  stop("return must be one of all, good or bad")
  #}

  singletons <- which(table(dge$samples$group) == 1)
  paired <- which(table(dge$samples$group) >= 2)

  # filter samples with technical replicates
  keep <- which(dge$samples$group %in% (names(paired)))
  dge.filter <- dge[,keep]

  corrs <- lapply(unique(names(paired)), function(x){
    data <- as.data.frame(dge$counts[,dge$samples$group==as.character(x)])
    fit <- stats::lm(data[,1] ~ data[,2])
    adj.r2 <- summary(fit)$adj.r.squared
    corr <- sqrt(adj.r2)
  })

  names(corrs) <- unique(names(paired))
  corrs <- unlist(corrs)

  # good corrs
  good.corrs <- corrs[which(corrs > corr.thresh)]
  # bad corrs
  bad.corrs <- corrs[which(corrs < corr.thresh)]

  if(return == "all"){
    return(corrs)
  }
  if(return == "bad"){
    return(bad.corrs)
  }
  if(return == "good"){
    return(good.corrs)
  }
}


