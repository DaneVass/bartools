# fixes issues with no visible binding for global variables. See here
# https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887
utils::globalVariables(c(".", "group", "percent", "value", "variable", "Var1", 
                         "Var2", "Sample", "Proportion", "Color", "barcode", 
                         "Group", "Barcodes", "Count", "rowid", "name", 
                         "test.dge", "BC.count", "model.matrix", 
                         "makeContrasts", "arrange", "FDR", "write.table", 
                         "logFC", "PValue", "category", "tail", "Freq", 
                         "num.barcodes","neglog10pval", "enriched", "cluster",
                         "Percentage", "n", "Barcode", "bc.umi.count", 
                         "cellid", "count", "count_max", "cumsum_count_max", 
                         "frac", "freq", "max_count", 
                         "number_of_lineage_barcodes", "proportion"))
