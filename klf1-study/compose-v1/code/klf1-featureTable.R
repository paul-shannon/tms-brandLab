library(trena)
library(ghdb)
library(igvR)
targetGene <- "KLF1"   # coded on the negative (reverse) strand
tss <- 12887201 #  chr19:12,887,201

ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg38")

tbl.fimo <- get(load("data/tbl.fimo.KLF1.chr19:12637201-13137201.RData"))
dim(tbl.fimo)    # 1576717       9
tbl.fimo$tssDistance <- tbl.fimo$end - tss  # upstream is positive
dim(tbl.fimo)

ft$setFundamentalRegions(tbl.fimo)
dim(ft$getTable())  # 157617 10

gh <- GeneHancerDB(host="ghdb", port=5432)
allTissues <- listTissues(gh, "KLF1")
tissues <-   c(grep("myeloid", allTissues, ignore.case=TRUE, v=TRUE),
               grep("eryth", allTissues, ignore.case=TRUE, v=TRUE),
               grep("blood", allTissues, ignore.case=TRUE, v=TRUE))
tbl.gh <- retrieveEnhancersFromDatabase(gh, "KLF1", tissues)
  # combined score has too much spread: 12 - 253.  compress
tbl.gh$combinedscore <- log2(tbl.gh$combinedscore)
dim(tbl.gh) # 3 16

feature.guide <- list(gh="combinedscore")
default.values <- list(gh=0)
ft$addRegionFeature(tbl.gh, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()
dim(tbl) # 1576717      11


mtx <- get(load("data/brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"))
dim(mtx) # 27171 28
fivenum(mtx["KLF1",])   # 7.89: 11.76

candidateTFs <- sort(unique(intersect(tbl.fimo$tf, rownames(mtx))))
length(candidateTFs) # 611
deleters <- grep("^KLF1$", candidateTFs)
if(length(deleters) > 0)
    candidateTFs <- candidateTFs[-deleters]
length(candidateTFs) # 610

solver <- EnsembleSolver(mtx,
                         targetGene=targetGene,
                         candidateRegulators=candidateTFs,
                         solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                         geneCutoff=1.0)
tbl.out <- suppressWarnings(run(solver))
new.order <- order(abs(tbl.out$spearmanCoeff), decreasing=TRUE)
tbl.out <- tbl.out[new.order,]
head(tbl.out)
dim(tbl.out)  # 602 7

ft$addGeneFeature(tbl.out[, c("gene", "rfScore")],
                  feature.name="rfScore", default.value=0)
ft$addGeneFeature(tbl.out[, c("gene", "spearmanCoeff")],
                  feature.name="spearmanCoeff", default.value=0)
ft$addGeneFeature(tbl.out[, c("gene", "xgboost")],
                  feature.name="xgboost", default.value=0)
ft$addGeneFeature(tbl.out[, c("gene", "betaLasso")],
                  feature.name="lasso", default.value=0)


tbl <- ft$getTable()

tbl.atac <- get(load("data/tbl.atacMerged.RData"))
dim(tbl.atac)
head(tbl.atac)

feature.guide <- list(atac="width")
default.values <- list(atac=0)
ft$addRegionFeature(tbl.atac, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()

tbl.gwas.raw <- get(load("data/tbl.gwas-vuckovic.RData"))
rownames(tbl.gwas.raw) <- NULL
coi <- c("CHR_ID", "CHR_POS", "PVALUE_MLOG", "cellType", "SNPS", "MAPPED_TRAIT")
tbl.gwas <- tbl.gwas.raw[, coi]
dim(tbl.gwas)  # 8403 6
head(tbl.gwas)
table(tbl.gwas$MAPPED_TRAIT)
dim(tbl.gwas)   # 8403 6
tbl.gwas <- subset(tbl.gwas, !is.na(CHR_POS))
dim(tbl.gwas)   # 8378 6

feature.guide <- list(gwasScore="PVALUE_MLOG", rsid="SNPS")
default.values <- list(gwasScore=0, rsid="")
tbl.gwas$start <- tbl.gwas$CHR_POS
tbl.gwas$end <- tbl.gwas$start + 1
tbl.gwas$chrom <- paste0("chr", tbl.gwas$CHR_ID)

ft$addRegionFeature(tbl.gwas, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()
igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "KLF1")
zoomOut(igv)
fivenum(tbl$rfScore)   # 0.000000e+00 4.893889e-05 1.892579e-03 3.329228e-02 9.294119e-01
fivenum(tbl$xgboost)   # 0.000000 0.000000 0.000000 0.000000 0.354464
fivenum(tbl$lasso)     # -0.003012231  0.000000000  0.000000000  0.000000000  0.303879632

tssThreshold <- 10000
tbl.strong.rf <- subset(tbl, atac > 0 & rfScore > 0.5 & abs(tssDistance) < tssThreshold)
dim(tbl.strong.rf)
sort(table(tbl.strong.rf$tf), decreasing=TRUE)
tbl.strong.lasso <- subset(tbl, atac > 0 & lasso > 0.1 & abs(tssDistance) < tssThreshold)
dim(tbl.strong.lasso)
sort(table(tbl.strong.lasso$tf), decreasing=TRUE)
tbl.strong.boost <- subset(tbl, atac > 0 & xgboost > 0.1 & abs(tssDistance) < tssThreshold)
dim(tbl.strong.boost)
sort(table(tbl.strong.boost$tf))

# having looked at the strongest candidates for each of these three methods,
# random forest, xgboost, and lasso, now combine them, and display some tracks
#
tbl.oi <- rbind(tbl.strong.rf, tbl.strong.lasso, tbl.strong.boost)
dim(tbl.oi)   # 121 18
tbl.oi <- unique(tbl.oi)
dim(tbl.oi)   # 68 18
new.order <- order(tbl.oi$start)
tbl.oi <- tbl.oi[new.order,]

sort(table(tbl.oi$tf), decreasing=TRUE)
  # E2F4 SREBF1  ARNT2   E2F7  GATA1  GFI1B
  #   51      7      5      2      2      1
track <- DataFrameQuantitativeTrack("", tbl.sub[, c("chrom", "start", "end", "combinedscore")],
                                    autoscale=TRUE)
displayTrack(igv, track)

tbl.strong <- subset(tbl, atac > 0 & abs(spearmanCoeff) > 0.8 & abs(tssDistance) < 2000)
dim(tbl.strong)
length(unique(tbl.strong$tf))
coi <- c()

tfs <- names(sort(table(tbl.sub.2$tf), decreasing=TRUE))
coi <- c("chrom", "start", "end", "spearmanCoeff", "tf")
tbl.oi <- tbl.sub.2[, ..coi]
