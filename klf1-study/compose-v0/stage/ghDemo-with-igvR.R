suppressWarnings(
   suppressMessages({
      library(ghdb)
      library(igvR)
      library(RPostgres)
   }))

gh <- GeneHancerDB(host="ghdb", port=5432)
allTissues <- listTissues(gh, "KLF1")
tissues <-   c(grep("myeloid", allTissues, ignore.case=TRUE, v=TRUE),
               grep("eryth", allTissues, ignore.case=TRUE, v=TRUE),
               grep("blood", allTissues, ignore.case=TRUE, v=TRUE))
tbl <- retrieveEnhancersFromDatabase(gh, "KLF1", tissues)
dim(tbl)

   #------------------------------------------------------------
   # now visualize the reported enhancers/promoters using
   # igvR.  there is one distal enhancer region which may
   # prove to be interesting
   # igvR opens the interactive genome view in a new window
   # or tab in your web browser.  RStudio continues to run
   # uninterrupted
   #------------------------------------------------------------

igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "KLF1")
zoomOut(igv)

   #------------------------------------------------------------
   # extract relevant columns from the genehancer table
   #------------------------------------------------------------

tbl.sub <- tbl[, c("chrom", "start", "end", "combinedscore")]

   #------------------------------------------------------------
   # the score seems topheavy.  reduce to more modest values
   #------------------------------------------------------------

tbl.sub$combinedscore  <- log2(tbl.sub$combinedscore)
track <- DataFrameQuantitativeTrack("gh", tbl.sub, autoscale=TRUE)
displayTrack(igv, track)



