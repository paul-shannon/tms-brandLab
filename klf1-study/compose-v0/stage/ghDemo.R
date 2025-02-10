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
View(tbl)

print(with(tbl, 1+end-start))

igv <- igvR()
setGenome(igv, "hg38")
tbl.sub <- tbl[, c("chrom", "start", "end", "combinedscore")]
tbl.sub$combinedscore <- log2(tbl.sub$combinedscore)
margin = 3000
klf1.region <- with(tbl.sub, sprintf("%s:%d-%d", chrom=chrom[1], min(start) - margin, max(end) + margin))
showGenomicRegion(igv, klf1.region)
track <- DataFrameQuantitativeTrack("gh", tbl.sub, autoscale=TRUE)
displayTrack(igv, track)






