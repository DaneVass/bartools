#' plotBarcodeRegression
#'
#' Generate a linear regression scatterplot for two sets of sample counts.
#'
#' @param dgeObject DGEList object with barcode counts.
#' @param sample1 Name of sample 1 (string).
#' @param sample2 Name of sample 2 (string).
#' @param title Optional, title of plot (string).
#' @param trendline Include linear trendline using `stat_smooth()` (boolean). Default = `TRUE`.
#' @param rug Include geom_rug density information on the axes (boolean). Defaule = `FALSE`.
#' @param trans Optional, the name of a transformation object or the object itself.
#'
#' @return Returns a scatterplot of the two samples
#' @export
#' @examples
#' data(test.dge)
#' plotBarcodeRegression(test.dge, sample1 = "T0-1", sample2 = "T0-2")

plotBarcodeRegression <-
  function(dgeObject,
           sample1 = NULL,
           sample2 = NULL,
           title = NULL,
           trendline = TRUE,
           trans = NULL,
           rug = FALSE) {
    inputChecks(dgeObject, samples = c(sample1, sample2))

    counts <- as.data.frame(dgeObject$counts)

    if (is.null(sample1) | is.null(sample2)) {
      stop("Must provide two samples to plot.")
    }

    plot.dat <- counts[, c(sample1, sample2)]

    # plot fit data
    # generate lm fit
    fit <- stats::lm(plot.dat[, 1] ~ plot.dat[, 2])

    p <-
      ggplot2::ggplot(fit$model, ggplot2::aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(title) +
      ggplot2::labs(y = colnames(plot.dat)[1], x = colnames(plot.dat)[2])

    if (!is.null(trans)) {
      p <- p + ggplot2::scale_x_continuous(trans = trans)
      p <- p + ggplot2::scale_y_continuous(trans = trans)
    }

    if (isTRUE(trendline)) {
      p <- p + ggplot2::stat_smooth(
        formula = y ~ x,
        method = "lm",
        col = "blue",
        se = F
      )
    }

    if (isTRUE(rug)) {
      p <-
        p + ggplot2::geom_rug(col = grDevices::rgb(.5, 0, 0, alpha = .2))
    }

    if (is.null(title)) {
      title <-
        paste("Regression plot:",
              as.character(sample1),
              "vs",
              as.character(sample2))
    }
    p <-
      p + ggplot2::ggtitle(
        title,
        subtitle = paste(
          "Adj R2 = ",
          signif(summary(fit)$adj.r.squared, 3),
          " ",
          " Intercept = ",
          signif(fit$coef[[1]], 3),
          " ",
          " Slope = ",
          signif(fit$coef[[2]], 3),
          " ",
          " P = ",
          summary(fit)$coef[2, 4]
        )
      )

    return(p)
  }
