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
#' @param counts DGElist list object with counts and samples dataframe
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



compareAbundance <- function(counts, meta, condition1,condition2,
                             pval.cutoff = 0.05, logFC.cutoff = 2,
                             filename = NULL) {

  cdt <- as.factor(counts$samples[[meta]])

  # barcode filtering using cpm > 0.1 in at least two samples
  cpm <- edgeR::cpm(counts$counts)
  countCheck <- cpm > 0.1
  keep <- which(rowSums(countCheck) >= 2)
  counts$counts <- counts$counts[keep,]

  dge = DGEList(counts=counts$counts)
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
  result <- as.data.frame(topTags(lrt, n = nrow(counts$counts)))
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
