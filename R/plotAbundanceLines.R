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
#' @param dge DGEList object containing raw counts of barcodes
#' @param condition sample condition of interest
#' @param conditionNames vector of size 2. Gives the name of the two conditions in conditions to plot
#' @param keep percentage of highest abundant barcode to keep. Rest of barcodes is filtered and not used in plots.
#' @param plotType which barcodes are plotted. 'DEG' plot the 10 most differentially abundant barcodes.
#' 'counts' plots the 5 top highest abundant barcodes for each condition.
#' 'log2FC' plots the barcodes with the highest absolute log2FC (number of barcodes to plot is given by nBarcodes)
#' @param nBarcodes number of barcodes to plot when using 'log2FC' plot type.
#' @param title desired plot title
#'
#' @return Returns a lineplot
#'
#' @export
#'
#' @examples
#' data(test.dge)
#' plotAbundanceLines(dge = test.dge, condition = test.dge$samples$group,
#' conditionNames = c("T0","10_High_dose"), plotType = 'counts', title = 'test')

plotAbundanceLines <-
  function(dge,
           condition,
           conditionNames,
           plotType = "DEG",
           keep = 0.9,
           nBarcodes = 10,
           title = "") {
    counts <- t(as.data.frame(dge$counts))

    # compute the mean for each barcode at each condition for the different mice
    df <-
      stats::aggregate(as.data.frame(counts), list(condition), stats::median)
    rownames(df) <- df$Group.1
    df$Group.1 <- NULL

    # keep the most expressed barcodes (%keep), remove rows with 0 values, then compute the log2 fold change. Only nBarcodes with the highest absolute FC are plotted.
    df <- as.data.frame(t(df))
    df <- df[names(df) %in% conditionNames]
    filtered_df <-
      utils::head(df[order(df[, 1], df[, 2], decreasing = T), ], nrow(df) * keep)

    ##Go through each row and determine if a value is zero
    row_sub = apply(filtered_df, 1, function(row)
      all(row != 0))
    ##Subset as usual
    filtered_df <- filtered_df[row_sub, ]

    if (plotType == "log2FC") {
      # compute the absolute log2 fold change for all barcodes for conditions
      filtered_df$abslog2FC <-
        abs(log2(filtered_df[, 1] / filtered_df[, 2]))
      filtered_df <-
        utils::head(filtered_df[order(filtered_df$abslog2FC, decreasing = T), ], nBarcodes)
      filtered_df$abslog2FC <- NULL
      legend_title <- "Highest log2FC barcodes"
    }

    if (plotType == "counts") {
      highest_barcodes <-
        c(rownames(utils::head(df[order(df[, 1], decreasing = T), ], 5)),
          rownames(utils::head(df[order(df[, 2], decreasing = T), ], 5)))
      filtered_df <- df[rownames(df) %in% highest_barcodes, ]
      legend_title <- "Highest count barcodes"
    }

    if (plotType == "DEG") {
      DEG_samples <-
        rownames(dge$samples[condition %in% conditionNames, ])
      DEG_group <- condition[condition %in% conditionNames]

      DEG_counts <- dge$counts[, colnames(dge$counts) %in% DEG_samples]
      DGE.obj <-
        edgeR::DGEList(counts = DEG_counts, group = DEG_group)

      # data filtering using cpm > 1
      cpm <- edgeR::cpm(DGE.obj)
      countCheck <- cpm > 1
      keep <- which(rowSums(countCheck) >= 2)
      DGE.obj.filtered <- DGE.obj[keep,]

      # normalisation - trimmed mean value
      DEG.obj.final <-
        edgeR::calcNormFactors(DGE.obj.filtered, method = "TMM")
      DEG.obj.final <-
        edgeR::estimateCommonDisp(DEG.obj.final, verbose = TRUE)
      DEG.obj.final <- edgeR::estimateTagwiseDisp(DEG.obj.final)
      DEG.res <- edgeR::exactTest(DEG.obj.final, pair = conditionNames)
      DEG.barcodes <-
        rownames(utils::head(DEG.res$table[order(DEG.res$table$PValue, decreasing = F), ], nBarcodes))
      filtered_df <- df[rownames(df) %in% DEG.barcodes, ]
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
