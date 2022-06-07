#' plotBarcodeRegression
#'
#' Generate a linear regression scatterplot for two sets of sample counts.
#'
#' @param dge DGEList object containing grouping variable to fit linear model
#' @param samp1 name of sample 1. must be one of colnames(dge)
#' @param samp2 name of sample 2. must be one of colnames(dge)
#' @param title desired name of output plot
#' @param trendline Logical. Include linear trendline using stat_smooth()
#' @param rug Logical. Include geom_rug density information on the axes?
#' @param trans From ggplot2. For continuous scales, the name of a transformation object or the object itself.
#'
#' @return Returns a scatterplot of the two samples
#' @export
#' @examples
#' data(test.dge)
#' plotBarcodeRegression(dge = test.dge, samp1 = "T0-1", samp2 = "T0-2")

plotBarcodeRegression <- function(dge, samp1 = NULL, samp2 = NULL, title = NULL, trendline = T, trans = NULL, rug = F) {
  dat <- as.data.frame(dge$counts)
  if( is.null(samp1) | is.null(samp2)) {
    print("must provide two samples to plot.")
    stop()
  }

  plot.dat <- as.data.frame(dat[,c(samp1,samp2)])

  # plot fit data
  # generate lm fit
  fit <- stats::lm(plot.dat[,1] ~ plot.dat[,2])

  p <- ggplot2::ggplot(fit$model, ggplot2::aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    ggplot2::labs(y = colnames(plot.dat)[1], x = colnames(plot.dat)[2])

  if(!is.null(trans)){
    p <- p + ggplot2::scale_x_continuous(trans = trans)
    p <- p + ggplot2::scale_y_continuous(trans = trans)
  }

  if(isTRUE(trendline)){
    p <- p + ggplot2::stat_smooth(method = "lm", col = "blue", se = F)
  }

  if(isTRUE(rug)){
    p <- p + ggplot2::geom_rug(col=grDevices::rgb(.5,0,0,alpha=.2))
  }

  if(is.null(title)){
    title <- paste("Regression plot:", as.character(samp1), "vs", as.character(samp2))
    p <- p + ggplot2::ggtitle(title, subtitle = paste("Adj R2 = ", signif(summary(fit)$adj.r.squared, 3), " ",
                                                      " Intercept = ", signif(fit$coef[[1]], 3), " ",
                                                      " Slope = ", signif(fit$coef[[2]], 3), " ",
                                                      " P = ", summary(fit)$coef[2,4]))
  } else {
    p <- p + ggplot2::ggtitle(title, subtitle = paste("Adj R2 = ", signif(summary(fit)$adj.r.squared, 3), " ",
                                                      " Intercept = ", signif(fit$coef[[1]], 3), " ",
                                                      " Slope = ", signif(fit$coef[[2]], 3), " ",
                                                      " P = ", summary(fit)$coef[2,4]))
  }
  return(p)
}
