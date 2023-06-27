#' compareAbundance
#'
#' @title
#' Barcode abundance comparison
#'
#' @description
#' Takes a DGElist list object with counts and samples dataframe.
#' Finds differentially abundant barcodes between two experimental conditions,
#' using the likelihood ratio test (LRT) based on the fitted generalized linear 
#' model (GLM) to compare two nested models, typically a reduced model and a 
#' full model, to determine if the additional terms in the full model 
#' significantly improve the model fit.  Additionally, it computes the p-values
#' associated with the test statistic for each barcode.
#' It save the results in a csv file and returns a volcano plot.
#'
#' @param dge.obj DGElist list object with counts and samples dataframe
#' @param meta Name of metadata field of interest for the comparison
#' @param condition1 Name of reference condition in meta
#' @param condition2 Name of second condition in meta
#' @param pval.cutoff Pvalue threshold for the volcano plot
#' @param logFC.cutoff logFC threshold for the volcano plot
#' @param filename Filename for the results of the likelihood ratio test
#' @return Returns a volcano plot and a csv file with differentially abundant 
#' barcode results
#' @export
#' @examples
#' data(test.dge)
#' compareAbundance(test.dge,"Treatment", "Vehicle","Low_dose",0.01,2)       



compareAbundance <- function(dge.obj, meta, condition1,condition2,
                             pval.cutoff = 0.05, logFC.cutoff = 2,
                             filename = NULL) {
  
  cdt <- as.factor(dge.obj$samples[[meta]])
  
  # barcode filtering using cpm > 0.1 in at least two samples
  cpm <- edgeR::cpm(dge.obj$counts)
  countCheck <- cpm > 0.1
  keep <- which(rowSums(countCheck) >= 2)
  dge.obj$counts <- dge.obj$counts[keep,]
  
  dge = DGEList(counts=dge.obj$counts)
  dge = calcNormFactors(dge, method='TMM') 
  design <- model.matrix(~cdt)
  colnames(design) <- levels(cdt)
  
  dge = estimateDisp(dge, design)
  
  # fits a generalized linear model (GLM) to the DGEList object
  fit <- glmFit(dge, design)

  # creates a contrast matrix to specify the comparison of interest
  pair_vector = sprintf("%s-%s", condition1, condition2) # Samples to be compared
  pair_contrast = makeContrasts(contrasts=pair_vector, levels=design) 
  
  # performs the likelihood ratio test based on the fitted model fit with the specified contrast
  lrt = glmLRT(fit, contrast=pair_contrast) # Likelihood ratio test
  
  #print(summary(decideTestsDGE(lrt)))

  # save results in a csv file ordered by FDR
  result <- as.data.frame(topTags(lrt, n = nrow(dge.obj$counts)))
  result <- dplyr::arrange(result, FDR)
  if (!is.null(filename)) {
    write.table(result, file = filename, sep="\t", row.names = T)
  }
  
  # Add a column to categorize the points based on the conditions
  result$category <- ifelse(result$PValue < pval.cutoff & result$logFC > logFC.cutoff, "Up",
                            ifelse(result$PValue < pval.cutoff & result$logFC < -logFC.cutoff, "Down", "Others"))
  result$label <- ifelse(result$PValue < pval.cutoff & abs(result$logFC) > logFC.cutoff, rownames(result), "")

  # Create a volcano plot
  plot <- ggplot(result, aes(x = logFC, y = -log10(PValue), color = category)) +
    geom_point(alpha = 0.6) +
    geom_text(aes(x = logFC, y = -log10(PValue)),label = result$label,vjust=1,parse = TRUE, show.legend = FALSE) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue"),
                       guide = guide_legend(title = NULL)) +
    labs(x = "logFC", y = "-log10(PValue)", title = paste(condition2,"vs", condition1,"differentially abundant barcodes", sep =" ")) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    scale_x_continuous(expand = c(0.1, 0.1))  # expand axes to see all BC labels
  
  # Return volcano plot of diffentially abundant barcodes
  return(plot)
}


# compareAbundance <- function(counts, condition, condition_names){
#   
#   DEG_samples <- rownames(counts$samples[condition %in% condition_names,])
#   DEG_group <- condition[condition %in% condition_names]
#   
#   DEG_counts <- counts$counts[,colnames(counts$counts) %in% DEG_samples]
#   
#   DGE.obj <- edgeR::DGEList(counts = DEG_counts, group =DEG_group)
#   
#   # data filtering using cpm > 1
#   cpm <- edgeR::cpm(DGE.obj)
#   countCheck <- cpm > 1
#   keep <- which(rowSums(countCheck) >= 2)
#   DGE.obj.filtered <- DGE.obj[keep,]
#   
#   # normalisation - trimmed mean value
#   DEG.obj.final <- edgeR::calcNormFactors(DGE.obj.filtered, method="TMM")
#   
#   DEG.obj.final <- edgeR::estimateCommonDisp(DEG.obj.final, verbose=TRUE)
#   DEG.obj.final <- edgeR::estimateTagwiseDisp(DEG.obj.final, trend="none")
#   
#   DEG.res <- edgeR::exactTest(DEG.obj.final, pair = condition_names)
#   DEG.res <- DEG.res$table[order(DEG.res$table$PValue,decreasing = F),]
#   return(DEG.res)
# }
