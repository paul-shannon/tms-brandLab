knitr::opts_chunk$set(echo = TRUE)

suppressWarnings(
  suppressMessages({
     library(trena)
     library(ghdb)
     library(igvR)
     }))

targetGene <- "KLF1"   # coded on the reverse strand
tss <- 12887201        #  chr19:12,887,201
ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg38")


#fimo.file <- "/home/rstudio/data/tbl.fimo.KLF1.chr19:12637201-13137201.RData"
fimo.file <- "./compose-v1/data/tbl.fimo.KLF1.chr19:12637201-13137201.RData"
tbl.fimo <- get(load(fimo.file))
dim(tbl.fimo)
tbl.fimo$tssDistance <- tbl.fimo$end - tss  # upstream is positive
dim(tbl.fimo)

ft$setFundamentalRegions(tbl.fimo)
dim(ft$getTable())

gh <- GeneHancerDB(host="ghdb", port=5432)
allTissues <- listTissues(gh, "KLF1")
tissues <-   c(grep("myeloid", allTissues, ignore.case=TRUE, v=TRUE),
               grep("eryth",   allTissues, ignore.case=TRUE, v=TRUE),
               grep("blood",   allTissues, ignore.case=TRUE, v=TRUE))
tbl.gh <- retrieveEnhancersFromDatabase(gh, "KLF1", tissues)
  # combined score has too much spread: 12 - 253.  compress
tbl.gh$combinedscore <- log2(tbl.gh$combinedscore)
dim(tbl.gh)

feature.guide <- list(gh="combinedscore")
default.values <- list(gh=0)
ft$addRegionFeature(tbl.gh, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()
dim(tbl)
colnames(tbl)

mtx <- get(load("/home/rstudio/data/brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"))
dim(mtx)

fivenum(mtx["KLF1",])

candidateTFs <- sort(unique(intersect(tbl.fimo$tf, rownames(mtx))))
length(candidateTFs)
deleters <- grep("^KLF1$", candidateTFs)
if(length(deleters) > 0)
    candidateTFs <- candidateTFs[-deleters]
length(candidateTFs)

solver <- EnsembleSolver(mtx,
                         targetGene=targetGene,
                         candidateRegulators=candidateTFs,
                         solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                         geneCutoff=1.0)
tbl.out <- suppressWarnings(run(solver))
new.order <- order(abs(tbl.out$spearmanCoeff), decreasing=TRUE)
tbl.out <- tbl.out[new.order,]
head(tbl.out)
dim(tbl.out)

ft$addGeneFeature(tbl.out[, c("gene", "rfScore")],
                  feature.name="rfScore", default.value=0)
ft$addGeneFeature(tbl.out[, c("gene", "spearmanCoeff")],
                  feature.name="spearmanCoeff", default.value=0)
ft$addGeneFeature(tbl.out[, c("gene", "xgboost")],
                  feature.name="xgboost", default.value=0)
ft$addGeneFeature(tbl.out[, c("gene", "betaLasso")],
                  feature.name="lasso", default.value=0)

tbl.atac <- get(load("/home/rstudio/data/tbl.atacMerged.RData"))
dim(tbl.atac)
head(tbl.atac)

feature.guide <- list(atac="width")
default.values <- list(atac=0)
ft$addRegionFeature(tbl.atac, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()
colnames(tbl)

tssThreshold <- 10000
tbl.strong.rf <- subset(tbl, atac > 0 & rfScore > 0.5 & abs(tssDistance) < tssThreshold & gh > 0)
dim(tbl.strong.rf)
sort(table(tbl.strong.rf$tf), decreasing=TRUE)

tbl.strong.lasso <- subset(tbl, atac > 0 & lasso > 0.1 & abs(tssDistance) < tssThreshold & gh > 0)
dim(tbl.strong.lasso)
sort(table(tbl.strong.lasso$tf), decreasing=TRUE)

tbl.strong.boost <- subset(tbl, atac > 0 & xgboost > 0.1 & abs(tssDistance) < tssThreshold & gh > 0)
dim(tbl.strong.boost)
sort(table(tbl.strong.boost$tf))

tbl.strong.spearman <- subset(tbl, atac > 0 & abs(spearmanCoeff) > 0.92 & gh > 0 &
                                abs(tssDistance) < tssThreshold)
sort(table(tbl.strong.spearman$tf))

tbl.c <- rbind(tbl.strong.rf, tbl.strong.lasso, tbl.strong.boost)
dim(tbl.c)
tbl.c <- unique(tbl.c)
dim(tbl.c)
new.order <- order(tbl.c$start)
tbl.c <- tbl.c[new.order,]
tf.tab <- sort(table(tbl.c$tf), decreasing=TRUE)
print(tf.tab)
tfs.c <- names(tf.tab)
tfs.spearman <- unique(tbl.strong.spearman$tf)

tfs.c

intersect(tfs.c, tfs.spearman)

setdiff(tfs.spearman, tfs.c)

igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "KLF1")
zoomOut(igv)
zoomOut(igv)
track <- DataFrameAnnotationTrack("GATA1", subset(tbl.c, tf=="GATA1"))
displayTrack(igv, track)
