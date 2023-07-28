## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bartools)

## ---- eval=FALSE--------------------------------------------------------------
#  pattern <- "([ACTG][ACTG][GC][AT][GC][ACTG][ACTG][AT][GC][AT]){3,6}"
#  fastq <- system.file("extdata", "test_extract_75bp_single-end.fastq.gz", package = "bartools", mustWork = T)
#  constant <- toupper("tgaccatgtacgattgacta")
#  test.extract <- bartools::extractbartoolseads(infile = fastq,
#                                barcode_pattern = pattern,
#                                constant = constant,
#                                yieldSize = 1e6)
#  test.extract

## ---- eval=FALSE--------------------------------------------------------------
#  test.map <- bartools::mapbartoolseads(reads = test.extract,
#                                        bowtie_index = "../data/bowtie/index",
#                                        mismatches = 1,
#                                        threads = 1,
#                                        prefix = "test_map")
#  test.map

## ---- eval=FALSE--------------------------------------------------------------
#  reference.fa <- system.file("extdata", "barcode_lib_reference_test.fasta", package = "bartools", mustWork = T)
#  test.map <- bartools::mapBarcodeReads(reads = test.extract,
#                                        bowtie_index = NULL,
#                                        reference_fasta = reference.fa,
#                                        mismatches = 1, # maximum 3 mismatches
#                                        threads = 2,
#                                        prefix = "test_map")
#  test.map

## ---- eval=FALSE--------------------------------------------------------------
#  test.counts <- utils::read.delim("./test_map_counts.csv", header = T, sep = ",", row.names = 1)
#  test.counts %>% tibble::rownames_to_column() %>% dplyr::arrange(dplyr::desc(mapped))

## -----------------------------------------------------------------------------
sessionInfo()

