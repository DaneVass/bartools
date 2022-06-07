#' processBarcodeLibrary
#'
#' process barcode reference file from raw sequencing datasets
#'
#' @param file path to starcode output for barcode library
#' @param samplename name of reference library. Will be prefixed to each barcode in the reference
#' @param cutoff rowsum cutoff defining barcodes to keep in the reference library
#' @param header Logical. Does the input file contain a header line
#' @param outdir Desired output directory to save library files
#'
#' @return returns data frame containing barcode and read count in reference
#'
#' @export
#'
#' @examples
#' load(system.file("extdata", "test_raw_lib.rda", package = "barista"))
#' processBarcodeLibrary(file = test.raw.lib, samplename = "Barcode", outdir = tempdir())

processBarcodeLibrary <- function(file = NULL, samplename = "Barcode", cutoff = 10, header = FALSE, outdir = tempdir()){
  # imports barcode reference file from starcode or similar data frame.
  # Barcode reference should have two columns as below
  # Barcode_ID Count
  # Barcode_1   100
  #
  if (is.null(file)){
    stop("Please include a valid path to starcode output")
  }

  if(class(file) == "data.frame"){
    barcodes <- file
  } else {
    barcodes <- data.table::fread(file)
  }

  if (header) {
    barcodes <- utils::tail(barcodes, -1) # get rid of the header line
  }

  # filter real barcodes
  colnames(barcodes) <- c("Barcode", "Raw_count")
  barcodes$Raw_count <- as.numeric(barcodes$Raw_count)
  barcodes.out <- barcodes[which(barcodes$Raw_count >= cutoff),] # set arbitrary cutoff of 10 as default

  # Generate rank and proportion information
  raw.total.count <- sum(barcodes$Raw_count)
  raw.proportion <- 100*(barcodes$Raw_count/raw.total.count)
  assertthat::are_equal(sum(raw.proportion), 100)


  filtered.total.count <- sum(barcodes.out$Raw_count)
  filtered.proportion <- 100*(barcodes.out$Raw_count/filtered.total.count)
  assertthat::are_equal(sum(filtered.proportion), 100)

  filtered.counts.starcode <- sum(barcodes[barcodes$Raw_count >= 100,2])
  proportion.starcode <- 100*(filtered.counts.starcode/raw.total.count)

  # output barcode library as fasta
  rank <- seq(1,length(rownames(barcodes.out)))
  rank <- paste(samplename, rank, sep = '')
  barcodes.out$Rank <- rank
  fasta.name <- paste(samplename,".fasta", sep = '')

  # output barcode library as fasta and txt
  message(paste("Writing library to", outdir))
  seqinr::write.fasta(sequences = as.list(barcodes.out$Barcode),
                      names = barcodes.out$Rank,
                      file.out = file.path(outdir,fasta.name), open = "w")
  utils::write.table(barcodes.out, file.path(outdir,paste(samplename, "_barcodes.txt", sep = "")), quote = F)

}


# get object name as string
getname <- function(v1) {
  deparse(substitute(v1))
}

# fasta output from dataframe function
writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,1], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,2]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

