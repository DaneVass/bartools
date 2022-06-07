#' Collapse technical replicates
#'
#' Collapse technical replicates in a DGEList object by mean or sum.
#' Modified from collapseReplicates function from DESeq2 package to accept DGEList objects and to allow collapsing replicated by mean or sum.
#'
#' @param object DGEList object containing raw or normalised barcode counts with replicate grouping info in object metadata
#' @param groupby a character vector of containing grouping information. Must be as long as the columns of object. Can pass a metadata column from object.
#' @param run optional, the names of each unique column in object. if provided, a new column runsCollapsed will be added to the colData which pastes together the names of run
#' @param renameCols boolean. whether to rename the columns of the returned object using the levels of the grouping factor
#' @param show_reps boolean. print replicate column names to console?
#' @param by collapse replicates by mean or sum
#'
#' @return Returns a dataframe of normalised counts for a sample
#' @export
#'
#' @examples
#' data(test.dge)
#' collapseReplicates(test.dge, groupby = test.dge$samples$group, by = "mean")
#'

collapseReplicates <- function (object, groupby, run, renameCols = TRUE, show_reps = TRUE, by = "mean") {

  # check inputs
  stopifnot(by != "mean" | by != "sum")

  if (!is.factor(groupby))
    groupby <- factor(groupby)
  groupby <- droplevels(groupby)
  stopifnot(length(groupby) == ncol(object))

  sp <- split(seq(along = groupby), groupby)
  if(isTRUE(show_reps)){
    print(sp)
  }

  # get obj class
  if(class(object)[1] == "DGEList"){
    countdata <- object$counts
  } else {
    stop("please supply a valid DGEList object as input")
  }

  # combine by mean or sum
  if(by == "mean"){
    countdata <- sapply(sp, function(i) rowMeans(countdata[,i, drop = FALSE]))
    mode(countdata) <- "integer"
  }
  if(by == "sum"){
    countdata <- sapply(sp, function(i) rowSums(countdata[,i, drop = FALSE]))
    mode(countdata) <- "integer"
  }

  # select columns to keep
  colsToKeep <- sapply(sp, `[`, 1)
  collapsed <- object[, colsToKeep]
  dimnames(countdata) <- dimnames(collapsed)
  collapsed$counts <- countdata

  if (!missing(run)) {
    stopifnot(length(groupby) == length(run))
    collapsed$samples$runsCollapsed <- sapply(sp, function(i) paste(run[i],
                                                                     collapse = ","))
  }

  # rename cols true by default
  if (renameCols) {
    colnames(collapsed) <- levels(groupby)
  }

  # output checks
  if(class(object)[1] == "DGEList"){
    stopifnot(sum(as.numeric(countdata)) == sum(as.numeric(collapsed$counts)))
  } else {
    stop("output is not a valid DGEList object")
  }
  collapsed
}


