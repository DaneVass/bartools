#' plotAbundanceLines
#'
#' Lineplot of the barcode abundances in two different experimental settings
#'
#' @title
#' Lineplot of barcode abundances in two conditions
#'
#' @description
#' Takes a dataframe of barcode counts,
#' computes the median abundance of each barcode for two specific conditions,
#' then do a line plot for both conditions
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param group Column name in sample metadata to group samples by (string).
#' @param conditions Names of 2 conditions in group column to compare (vector of strings).
#' @param keep percentage of highest abundant barcode to keep (decimal). Rest of barcodes is filtered and not used in plots. Default = `0.9`.
#' @param plotType (string) `DEG` plots the 10 most differentially abundant barcodes.
#' `counts` plots the 5 top highest abundant barcodes for each condition.
#' `log2FC` plots the barcodes with the highest absolute log2FC (number of barcodes to plot is given by nBarcodes).
#' Default = `DEG`.
#' @param nBarcodes Number of barcodes to plot when using 'log2FC' plot type (integer). Default = `10`.
#' @param title Optional, plot title (string).
#'
#' @return Returns a lineplot
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' plotAbundanceLines(test.dge, group = "group",
#' conditions = c("T0","10_High_dose"), plotType = "counts")

plotAbundanceLines <-
  function(dgeObject,
           group,
           conditions,
           plotType = "DEG",
           keep = 0.9,
           nBarcodes = 10,
           title = "") {

    inputChecks(dgeObject, conditions = conditions, groups = group)

    counts <- t(as.data.frame(dgeObject$counts))

    # compute the mean for each barcode at each condition for the different mice
    df <-
      stats::aggregate(as.data.frame(counts), list(group), stats::median)
    rownames(df) <- df$Group.1
    df$Group.1 <- NULL

    # keep the most expressed barcodes (%keep), remove rows with 0 values, then compute the log2 fold change. Only nBarcodes with the highest absolute FC are plotted.
    df <- as.data.frame(t(df))
    df <- df[names(df) %in% conditions]
    filtered_df <-
      utils::head(df[order(df[, 1], df[, 2], decreasing = T),], nrow(df) * keep)

    ##Go through each row and determine if a value is zero
    row_sub = apply(filtered_df, 1, function(row)
      all(row != 0))
    ##Subset as usual
    filtered_df <- filtered_df[row_sub,]

    if (plotType == "log2FC") {
      # compute the absolute log2 fold change for all barcodes for conditions
      filtered_df$abslog2FC <-
        abs(log2(filtered_df[, 1] / filtered_df[, 2]))
      filtered_df <-
        utils::head(filtered_df[order(filtered_df$abslog2FC, decreasing = T),], nBarcodes)
      filtered_df$abslog2FC <- NULL
      legend_title <- "Highest log2FC barcodes"
    }

    if (plotType == "counts") {
      highest_barcodes <-
        c(rownames(utils::head(df[order(df[, 1], decreasing = T),], 5)),
          rownames(utils::head(df[order(df[, 2], decreasing = T),], 5)))
      filtered_df <- df[rownames(df) %in% highest_barcodes,]
      legend_title <- "Highest count barcodes"
    }

    if (plotType == "DEG") {
      DEG_samples <-
        rownames(dgeObject$samples[group %in% conditions,])
      DEG_group <- group[group %in% conditions]

      DEG_counts <-
        dgeObject$counts[, colnames(dgeObject$counts) %in% DEG_samples]
      DGE.obj <-
        edgeR::DGEList(counts = DEG_counts, group = DEG_group)

      # data filtering using cpm > 1
      cpm <- edgeR::cpm(DGE.obj)
      countCheck <- cpm > 1
      keep <- which(rowSums(countCheck) >= 2)
      DGE.obj.filtered <- DGE.obj[keep, ]

      # normalisation - trimmed mean value
      DEG.obj.final <-
        edgeR::calcNormFactors(DGE.obj.filtered, method = "TMM")
      DEG.obj.final <-
        edgeR::estimateCommonDisp(DEG.obj.final, verbose = TRUE)
      DEG.obj.final <- edgeR::estimateTagwiseDisp(DEG.obj.final)
      DEG.res <-
        edgeR::exactTest(DEG.obj.final, pair = conditions)
      DEG.barcodes <-
        rownames(utils::head(DEG.res$table[order(DEG.res$table$PValue, decreasing = F),], nBarcodes))
      filtered_df <- df[rownames(df) %in% DEG.barcodes,]
      legend_title <- "Differentially abundant barcodes"
    }

    # plot data
    melted_df <- reshape2::melt(filtered_df, id.vars = NULL)
    melted_df$rowid <- rownames(filtered_df)
    ggplot2::ggplot(melted_df, ggplot2::aes(variable, value, group = factor(rowid))) +
      ggplot2::geom_line(ggplot2::aes(color = factor(rowid))) +
      ggplot2::xlab("") +
      ggplot2::ylab("Counts") +
      ggplot2::scale_colour_discrete(legend_title) +
      ggplot2::ggtitle(title) +
      ggplot2::theme_bw()
  }
