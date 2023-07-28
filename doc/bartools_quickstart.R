## ----setup--------------------------------------------------------------------
suppressWarnings(library(bartools))
knitr::opts_chunk$set(dev="png")

## -----------------------------------------------------------------------------
data(test.dge)

## ----eval = FALSE-------------------------------------------------------------
#  samplesheet <-
#    read.csv(
#      system.file(
#        "extdata",
#        "test_sampletable.csv",
#        package = "bartools",
#        mustWork = T
#      ),
#      header = T,
#      stringsAsFactors = F
#    )
#  samplesheet

## ----eval = FALSE-------------------------------------------------------------
#  dge <-
#    edgeR::readDGE(
#      files = samplesheet,
#      group = samplesheet$treatment,
#      labels = samplesheet$sample,
#      header = T
#    )

## -----------------------------------------------------------------------------
data(test.dge)
test.dge

## -----------------------------------------------------------------------------
# Remove rows with no data
thresholdCounts(test.dge, type = "absolute", threshold = 1, min.samps = 1, plot = T, group = "Treatment")
thresholdCounts(test.dge, type = "absolute", threshold = 10, min.samps = 1, plot = T, group = "Treatment")
thresholdCounts(test.dge, type = "absolute", threshold = 10, min.samps = 3, plot = T, group = "Treatment")

## -----------------------------------------------------------------------------
# Remove rows with no data
thresholdCounts(test.dge, type = "relative", threshold = 1e-10, min.samps = 1, plot = T, group = "Treatment")
thresholdCounts(test.dge, type = "relative", threshold = 1e-5, min.samps = 1, plot = T, group = "Treatment")
thresholdCounts(test.dge, type = "relative", threshold = 1e-5, min.samps = 3, plot = T, group = "Treatment")

## -----------------------------------------------------------------------------
dge.filtered <- thresholdCounts(test.dge, type = "absolute", threshold = 10, min.samps = 2, plot = F)

## -----------------------------------------------------------------------------
dge.cpmnorm <- normaliseCounts(dge.filtered, method = "CPM")

## -----------------------------------------------------------------------------
# raw counts per sample
plotReadCounts(dge.filtered$counts, group = dge.filtered$samples$Treatment)

## -----------------------------------------------------------------------------
# normalised counts per sample
plotReadCounts(dge.cpmnorm, group = dge.filtered$samples$Treatment)

## -----------------------------------------------------------------------------
# plot detected barcodes ordered by frequency in reference library
plotBarcodeCounts(dge.cpmnorm, log10 = F)
# plot log10 barcode counts
plotBarcodeCounts(dge.cpmnorm, log10 = T)
# order barcodes by count across samples
plotBarcodeCounts(dge.cpmnorm, log10 = T, order = T)

## ----message=FALSE------------------------------------------------------------
samps <- unique(dge.filtered$samples$group)
# only plot subset of samples
lapply(samps[4:6], function(x) {
  df <- dge.filtered[, dge.filtered$samples$group %in% as.character(x)]
  plotBarcodeRegression(
    df,
    samp1 = colnames(df)[[1]],
    samp2 = colnames(df)[[2]],
    rug = T,
    trans = "log1p"
  )
})

## -----------------------------------------------------------------------------
corrs <- calcReplicateCorr(dge.filtered, group = "Group")
corrs <- as.data.frame(corrs)
corrs

## -----------------------------------------------------------------------------
which(corrs < 0.9)

## -----------------------------------------------------------------------------
p <- ggplot(reshape2::melt(corrs, value.name = "Correlation"),
            aes(x = variable, y = Correlation)) +
  geom_boxplot(width = 0.4, outlier.size = 0) +
  geom_jitter(width = 0.02, ) +
  scale_y_continuous(trans = "log1p") +
  theme_classic()
p

plotBarcodeCorrelation(dge.filtered$counts, clustered = T, upper = T, method = "pearson")
plotBarcodeCorrelation(dge.filtered$counts, clustered = T, upper = T, method = "spearman")

## -----------------------------------------------------------------------------
dim(dge.filtered)

## -----------------------------------------------------------------------------
dge.filtered.collapsed <- collapseReplicates(
  dge.filtered,
  groupby = dge.filtered$samples$group,
  by = "mean",
  show_reps = F
)
dge.filtered.collapsed <- dge.filtered.collapsed[, rev(order(dge.filtered.collapsed$samples$Treatment))]

## -----------------------------------------------------------------------------
dim(dge.filtered.collapsed)

## -----------------------------------------------------------------------------
head(dge.filtered.collapsed)

## -----------------------------------------------------------------------------
plotBarcodeBubble(dge.filtered.collapsed$counts, proportion.cutoff = 10, labels = T)

## -----------------------------------------------------------------------------
plotOrderedBubble(counts.obj = dge.filtered.collapsed$counts, proportion.cutoff = 10, labels = T, orderSample = "T0", colorDominant = F, filterLow = T, samples = dge.filtered.collapsed$samples, group = "Treatment")

## -----------------------------------------------------------------------------
plotOrderedBubble(
  dge.filtered.collapsed$counts,
  proportion.cutoff = 10,
  labels = T,
  orderSample = "T0",
  colorDominant = T
)

## -----------------------------------------------------------------------------
plotBarcodeHistogram(dge.filtered.collapsed$counts,
                     sample = dge.filtered.collapsed$samples$group[[10]],
                     top = 50)

## -----------------------------------------------------------------------------
plotBarcodeTimeseries(dge.filtered.collapsed, top = 5)

## -----------------------------------------------------------------------------
plotBarcodePCA(dge.filtered.collapsed, intgroup = "Treatment")
plotBarcodePCA(
  dge.filtered.collapsed[, dge.filtered.collapsed$samples$Treatment %in% c("Baseline", "Vehicle", "High_dose")], intgroup = "Treatment"
)

## -----------------------------------------------------------------------------
plotBarcodeHeatmap(
  counts = cpm(dge.filtered.collapsed$counts),
  N = 5,
  show_bc = T,
  samples = dge.filtered.collapsed$samples,
  group = "Treatment"
)

## -----------------------------------------------------------------------------
plotBarcodeCumSum(dge.filtered.collapsed$counts, sample1 = "T0", samples = colnames(dge.filtered.collapsed$counts)[1:5])

## -----------------------------------------------------------------------------
top.bc <- getDominantBarcodes(dge.filtered.collapsed, pct.thresh = 5)
top.bc[1:5]

## -----------------------------------------------------------------------------
plotBarcodeBoxplot(dge.filtered.collapsed, barcodes = top.bc$T0, condition = c("T0", "Low_dose", "High_dose"))

## -----------------------------------------------------------------------------
top_barcodes <- calcPercentileBarcodes(dge.filtered.collapsed, percentile = 0.95)

top_barcodes$NumBarcodes
top_barcodes$TopBarcodeCounts$`6_High_dose`
top_barcodes$TopBarcodes$`6_High_dose`

## -----------------------------------------------------------------------------
plotDetectedBarcodes(
  dge.filtered.collapsed,
  percentile = 1,
  plot = T,
  group = "Treatment"
)
plotDetectedBarcodes(
  dge.filtered.collapsed,
  percentile = 0.95,
  plot = T,
  group = "Treatment"
)

## -----------------------------------------------------------------------------
diversity <- calcDivIndexes(dge.filtered.collapsed$counts)
diversity

## -----------------------------------------------------------------------------
qplot(diversity$name, diversity$shannon) + theme_bw() + coord_flip()

## -----------------------------------------------------------------------------
compareAbundance(dge.obj = dge.filtered.collapsed,
                 meta = "Treatment", 
                 condition1 = "Vehicle",
                 condition2 = "Low_dose",
                 pval.cutoff = 0.001,
                 logFC.cutoff = 5)


## -----------------------------------------------------------------------------
plotAbundanceLines(
  dge.filtered.collapsed,
  condition = dge.filtered.collapsed$samples$Treatment,
  condition_names = c("Vehicle", "High_dose"),
  plot_type = 'counts'
)

## -----------------------------------------------------------------------------
sessionInfo()

