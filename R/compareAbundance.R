#' compareAbundance
#'
#' Compare barcode abundance between 2 conditions using EdgeR
#'
#' @title
#' Barcode abundance comparison
#'
#' @description
#' Takes a dataframe of barcode counts,
#' computes the mean abundance of each barcode for two specific conditions,
#' then do a line plot for both conditions
#'
#' @param counts dataframe containing raw counts of barcodes
#' @param condition sample condition of interest
#' @param condition_names names of possible conditions in experiment
#' @return Returns a text file
#' @export
#' @examples
#' data(test.dge)
#' compareAbundance(test.dge, test.dge$samples$group, c("T0","T0b"))


compareAbundance <- function(counts, condition, condition_names){

  DEG_samples <- rownames(counts$samples[condition %in% condition_names,])
  DEG_group <- condition[condition %in% condition_names]

  DEG_counts <- counts$counts[,colnames(counts$counts) %in% DEG_samples]

  DGE.obj <- edgeR::DGEList(counts = DEG_counts, group =DEG_group)

  # data filtering using cpm > 1
  cpm <- edgeR::cpm(DGE.obj)
  countCheck <- cpm > 1
  keep <- which(rowSums(countCheck) >= 2)
  DGE.obj.filtered <- DGE.obj[keep,]

  # normalisation - trimmed mean value
  DEG.obj.final <- edgeR::calcNormFactors(DGE.obj.filtered, method="TMM")

  DEG.obj.final <- edgeR::estimateCommonDisp(DEG.obj.final, verbose=TRUE)
  DEG.obj.final <- edgeR::estimateTagwiseDisp(DEG.obj.final, trend="none")

  DEG.res <- edgeR::exactTest(DEG.obj.final, pair = condition_names)
  DEG.res <- DEG.res$table[order(DEG.res$table$PValue,decreasing = F),]
  return(DEG.res)
}
