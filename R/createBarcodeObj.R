#' createBarcodeObj
#'
#' generate a DGEList object containing barcode counts and sample metadata from raw counts files
#'
#' @param counts.dir directory containing counts files per sample
#' @param sampleanno annotation file detailing sample ID and groups
#'
#' @return A data frame containing barcode counts in rows per sample in columns.
#'
#' @import edgeR
#'

createBarcodeObj <- function(counts.dir = "counts/", sampleanno = NULL){
  counts <- edgeR::readDGE(files = list.files(counts.dir), path = counts.dir, group = sampleanno$group, labels = sampleanno$ID)
  counts$samples$replicate <- sampleanno$Replicate
  dim(counts)

  # raw.counts <- as.data.frame(counts$counts)
  #
  # # assert that samples are ordered correctly
  # for (i in rownames(counts$samples)){
  #   print(i)
  #   sample = paste("counts/",i,".counts.txt", sep='')
  #   index <- as.numeric(which(files == sample))
  #   df <- read.delim(files[index])
  #   if(sum(df[2]) == counts$samples[i,3]){
  #     print("CORRECT")
  #   } else {
  #     print("INCORRECT")
  #   }
  #
  # }
}

