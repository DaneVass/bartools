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
#' @param dgeObject DGEList object with barcode counts.
#' @param group Column name in sample metadata to group samples by (string).
#' @param condition1 Name of reference condition in group column (string).
#' @param condition2 Name of second condition in group column (string).
#' @param pval Pvalue threshold for the volcano plot (decimal). Default = `0.05`.
#' @param logFC logFC threshold for the volcano plot (decimal). Default = `2`.
#' @param filename Optional, filename for the results of the likelihood ratio test (string).
#'
#' @return Returns a volcano plot and a csv file with differentially abundant
#' barcode results
#' @export
#' @examples
#' data(test.dge)
#' compareAbundance(test.dge, group = "Treatment", condition1 = "Vehicle", condition2 = "Low_dose")

compareAbundance <-
  function(dgeObject,
           group,
           condition1,
           condition2,
           pval = 0.05,
           logFC = 2,
           filename = NULL) {
    inputChecks(dgeObject,
                groups = group,
                conditions = c(condition1, condition2))

    cdt <- as.factor(dgeObject$samples[[group]])

    # barcode filtering using cpm > 0.1 in at least two samples
    cpm <- edgeR::cpm(dgeObject$counts)
    countCheck <- cpm > 0.1
    keep <- which(rowSums(countCheck) >= 2)
    dgeObject$counts <- dgeObject$counts[keep,]

    dge = DGEList(counts = dgeObject$counts)
    dge = calcNormFactors(dge, method = 'TMM')
    design <- model.matrix(~ cdt)
    colnames(design) <- levels(cdt)

    dge = estimateDisp(dge, design)

    # fits a generalized linear model (GLM) to the DGEList object
    fit <- glmFit(dge, design)

    # creates a contrast matrix to specify the comparison of interest
    # conditions must be valid R variables, i.e. start with a character
    colnames(design) <- sapply(colnames(design), FUN = function(x) if (grepl("^[[:digit:]]", x)) paste0("char_", x) else x)
    condition1_variable <- if (grepl("^[[:digit:]]", condition1)) paste0("char_", condition1) else condition1
    condition2_variable <- if (grepl("^[[:digit:]]", condition2)) paste0("char_", condition2) else condition2

    pair_vector = sprintf("%s-%s", condition1_variable, condition2_variable) # Samples to be compared
    pair_contrast = makeContrasts(contrasts = pair_vector, levels = design)
    # performs the likelihood ratio test based on the fitted model fit with the specified contrast
    lrt = glmLRT(fit, contrast = pair_contrast) # Likelihood ratio test

    #print(summary(decideTestsDGE(lrt)))

    # save results in a csv file ordered by FDR
    result <-
      as.data.frame(topTags(lrt, n = nrow(dgeObject$counts)))
    result <- dplyr::arrange(result, FDR)
    if (!is.null(filename)) {
      write.table(result,
                  file = filename,
                  sep = "\t",
                  row.names = T)
    }

    # Add a column to categorize the points based on the conditions
    result$category <-
      ifelse(
        result$PValue < pval & result$logFC > logFC,
        "Up",
        ifelse(
          result$PValue < pval &
            result$logFC < -logFC,
          "Down",
          "Others"
        )
      )
    result$label <-
      ifelse(result$PValue < pval &
               abs(result$logFC) > logFC,
             rownames(result),
             "")

    # Create a volcano plot
    plot <-
      ggplot(result, aes(
        x = logFC,
        y = -log10(PValue),
        color = category
      )) +
      geom_point(alpha = 0.6) +
      geom_text(
        aes(x = logFC, y = -log10(PValue)),
        label = result$label,
        vjust = 1,
        parse = TRUE,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Others" = "gray"),
                         guide = guide_legend(title = NULL)) +
      labs(
        x = "logFC",
        y = "-log10(PValue)",
        title = paste(
          condition2,
          "vs",
          condition1,
          "differentially abundant barcodes",
          sep = " "
        )
      ) +
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
