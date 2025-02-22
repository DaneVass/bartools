# bartools 1.1.0

- Fix plotDetectedBarcodes and calcPercentileBarcodes to be inclusive, i.e. at least 1 barcode makes up top x percentile of a sample
- improve coloring of barcodes in bubble plots
- update plotDetectedBarcodes
- collapseReplicates and thresholdCounts update the lib.size column to the sum of the merged replicates or filtered sample.
- fix bug where small dots were plotted for barcodes that were not present in sample in plotOrderedBubble and plotBarcodeBubble

# bartools 1.0.0

- Updated function parameter inputs across whole package
- Updated quickstart vignette, added single-cell vignette
- Added functions aggregateBarcodes, filterBarcodes, plotUmiFilterThresholds, plotUmiPerBarcode, readBartabCounts, plotBarcodesPerCell for single-cell data and plotSampleCumSum for bulk data
- removed deprecated plotLibraryCumSum, plotLibraryDiversity, proportionalBubbleplot

# bartools 0.2.5

-   Added single cell QC plotting functions plotCellsPerGroup, plotMetrics, plotClusterEnrichment, plotCellsinClusters

# bartools 0.2.4

-   Updated plotting functions for plotBarcodeBoxplot, plotAbundanceLines

# bartools 0.2.3

-   Removed deprecated single cell workflows from package (workflow now in BARtab)
-   Updated plotting functions for plotOrderedBubble, plotBarcodeBubble, plotBarcodeHeatmap

# bartools 0.2.2

-   Added dose escalation test dataset to package 

# bartools 0.2.1

-   Updated plotting functions

# bartools 0.2.0

-   Updated plotDetectedBarcodes to color and reorder factors according to user input.
-   Updated quickstart vignette detailing a complete analysis workflow.

# bartools 0.1.0

-   Updated plotBarcodeTimeseries to take top n barcodes. Speeds up plotting by only focusing on barcodes of interest

# bartools 0.1.0

-   Included current versions of all plotting functions to package.
-   Added a `NEWS.md` file to track changes to the package.
