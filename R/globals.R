# fixes issues with no visible binding for global variables. See here
# https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887
utils::globalVariables(c(".", "group", "percent", "value", "variable", "Var1", "Var2", "Sample", 
                         "Proportion", "Color", "barcode", "Group", "Barcodes", "Count", "rowid",
                         "name", "test.dge", "BC.count", "model.matrix", "makeContrasts", "arrange", "FDR", "write.table",
                         "logFC", "PValue", "category"))

colors.npg <- c("#E64B35","#D85543","#CA5F52","#BC6960","#AE736F","#A07D7D","#92888C","#84929A","#769CA9","#68A6B7"
  ,"#5AB0C6","#4CBBD5","#45B8CD","#3EB6C6","#37B3BF","#30B1B8","#29AEB1","#22ACAA","#1BA9A3","#14A79C"
  ,"#0DA495","#06A28E","#009F87","#059987","#0A9287","#108B87","#158487","#1B7D87","#207687","#266F87"
  ,"#2B6887","#316187","#365A87","#3C5488","#4C5A87","#5D6086","#6D6785","#7E6D84","#8F7483","#9F7A83"
  ,"#B08182","#C18781","#D18E80","#E2947F","#F29B7F","#E89A83","#DE9988","#D4988D","#CA9792","#C09697"
  ,"#B6959B","#AC94A0","#A293A5","#9892AA","#8E91AF","#8491B4","#8596B5","#869CB6","#87A2B7","#88A8B9"
  ,"#89AEBA","#8BB3BB","#8CB9BC","#8DBFBE","#8EC5BF","#8FCBC0","#91D0C1","#97BDB0","#9EAA9E","#A5978D"
  ,"#AC847B","#B37169","#B95E58","#C04B46","#C73834","#CE2523","#D51211","#DB0000","#D30806","#CA110D"
  ,"#C21A13","#B9231A","#B12C20","#A83427","#A03D2D","#974634","#8F4F3A","#865841","#7E6148","#82664D"
  ,"#876B53","#8B7158","#90765E","#947B63","#998169","#9D866E","#A28B74","#A69179","#AB967F","#B09C85")
