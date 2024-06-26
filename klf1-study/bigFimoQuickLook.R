library(igvR)
igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "KLF1")
tbl.bf <- get(load("tbl.fimo.KLF1.chr19:12637201-13137201.RData"))
fivenum(tbl.bf$score)
p.max <- 1e-9
dim(subset(tbl.bf, p.value < p.max))
tbl.sub <- subset(tbl.bf, p.value < p.max)
track <- DataFrameAnnotationTrack("bf 1e-9", tbl.sub[, c("chrom", "start", "end")])
displayTrack(igv, track)

