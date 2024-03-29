% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotClusterEnrichment.R
\name{plotClusterEnrichment}
\alias{plotClusterEnrichment}
\title{Hypergeometric test for cluster enrichment, e.g. to test whether a barcode is enriched in any cluster.}
\usage{
plotClusterEnrichment(
  sc.obj = NULL,
  group = "barcode",
  factor = NULL,
  clusters = "seurat_clusters",
  threshold = 0.01,
  order = TRUE,
  plot = TRUE
)
}
\arguments{
\item{sc.obj}{Single cell Seurat or SingleCellExperiment object containing clusters and group metadata.}

\item{group}{a column of metadata (string). Default = `barcode`.}

\item{factor}{a level of group to test for enrichment per cluster, e.g. a specific barcode (string).}

\item{clusters}{a column of metadata defining the cluster identities of each cell (string). Default = `seurat_clusters`}

\item{threshold}{P-value threshold for hypergeometric test (decimal). Default = `0.01`.}

\item{order}{Order clusters in plot by -log10 P-value (boolean). Default = `TRUE`.}

\item{plot}{Plot hypergeometric test results or return plot data (boolean). Default = `TRUE`.}
}
\value{
Returns a histogram of hypergeometric test results or underlying plot data
}
\description{
Takes a single cell object, a grouping variable, a factor within the group to test,
and an ident class (i.e. clusters). Per level of ident performs hypergeometric testing
for enrichment of group factor
}
\details{
plotClusterEnrichment

Hypergeometric test for enrichment of specific classes of cells in
single cell RNA-seq clusters
}
