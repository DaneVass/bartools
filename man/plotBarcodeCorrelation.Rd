% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotBarcodeCorrelation.R
\name{plotBarcodeCorrelation}
\alias{plotBarcodeCorrelation}
\title{Sample to sample correlation}
\usage{
plotBarcodeCorrelation(
  dgeObject,
  method = "pearson",
  upper = TRUE,
  clustered = TRUE,
  title = "Sample Correlation"
)
}
\arguments{
\item{dgeObject}{DGEList object with barcode counts.}

\item{method}{Correlation metric to use (string). One of "pearson", "kendall", or "spearman". Refer to stats:cor for details. Default = `pearson`.}

\item{upper}{Plot only the upper half of the matrix (boolean). Default = `TRUE`.}

\item{clustered}{Cluster rows and columns of matrix (boolean). Default = `TRUE`.}

\item{title}{Title of plot (string). Default = `Sample Correlation - method`.}
}
\value{
Returns a heatmap of sample distances using desired clustering
}
\description{
Plot sample correlation between barcode sets / samples
}
\examples{
data(test.dge)
plotBarcodeCorrelation(test.dge)
}
