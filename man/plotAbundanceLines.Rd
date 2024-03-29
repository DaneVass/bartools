% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAbundanceLines.R
\name{plotAbundanceLines}
\alias{plotAbundanceLines}
\title{Lineplot of barcode abundances in two conditions}
\usage{
plotAbundanceLines(
  dgeObject,
  group,
  conditions,
  plotType = "DEG",
  keep = 0.9,
  nBarcodes = 10,
  title = ""
)
}
\arguments{
\item{dgeObject}{DGEList object with barcode counts.}

\item{group}{Column name in sample metadata to group samples by (string).}

\item{conditions}{Names of 2 conditions in group column to compare (vector of strings).}

\item{plotType}{(string) `DEG` plots the 10 most differentially abundant barcodes.
`counts` plots the 5 top highest abundant barcodes for each condition.
`log2FC` plots the barcodes with the highest absolute log2FC (number of barcodes to plot is given by nBarcodes).
Default = `DEG`.}

\item{keep}{percentage of highest abundant barcode to keep (decimal). Rest of barcodes is filtered and not used in plots. Default = `0.9`.}

\item{nBarcodes}{Number of barcodes to plot when using 'log2FC' plot type (integer). Default = `10`.}

\item{title}{Optional, plot title (string).}
}
\value{
Returns a lineplot
}
\description{
Takes a dataframe of barcode counts,
computes the median abundance of each barcode for two specific conditions,
then do a line plot for both conditions
}
\details{
plotAbundanceLines

Lineplot of the barcode abundances in two different experimental settings
}
\examples{
data(test.dge)
plotAbundanceLines(test.dge, group = "group",
conditions = c("T0","10_High_dose"), plotType = "counts")
}
