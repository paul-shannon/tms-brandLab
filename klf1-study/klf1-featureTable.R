library(trena)
library(ghdb)
targetGene <- "KLF1"   # coded on the negative (reverse) strand
tss <- 12887201 #  chr19:12,887,201

ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg38")

  # tbl.fimo <- get(load("../fimo/KLF1-20k-em5/fimo-5.RData"))
  # dim(tbl.fimo)  #   40664     9
tbl.fimo <- get(load("tbl.fimo.KLF1.chr19:12637201-13137201.RData"))
dim(tbl.fimo)    # 1576717       9
tbl.fimo$tssDistance <- tbl.fimo$end - tss  # upstream is positive
dim(tbl.fimo)

ft$setFundamentalRegions(tbl.fimo)
dim(ft$getTable())


gh <- GeneHancerDB()
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



mtx <- get(load("../featureTable/brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"))
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
tbl <- ft$getTable()

tbl.atac <- get(load("../featureTable/tbl.atacMerged.RData"))
dim(tbl.atac)
head(tbl.atac)

feature.guide <- list(atac="width")
default.values <- list(atac=0)
ft$addRegionFeature(tbl.atac, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()

tbl.gwas.raw <- get(load("../featureTable/tbl.gwas-vuckovic.RData"))
rownames(tbl.gwas.raw) <- NULL
coi <- c("CHR_ID", "CHR_POS", "PVALUE_MLOG", "cellType", "SNPS", "MAPPED_TRAIT")
tbl.gwas <- tbl.gwas.raw[, coi]
dim(tbl.gwas)  # 8403 5
head(tbl.gwas)
table(tbl.gwas$MAPPED_TRAIT)
dim(tbl.gwas)   # 8403 9
tbl.gwas <- subset(tbl.gwas, !is.na(CHR_POS))
dim(tbl.gwas)   # 8378 9

feature.guide <- list(gwasScore="PVALUE_MLOG", rsid="SNPS")
default.values <- list(gwasScore=0, rsid="")
tbl.gwas$start <- tbl.gwas$CHR_POS
tbl.gwas$end <- tbl.gwas$start + 1
tbl.gwas$chrom <- paste0("chr", tbl.gwas$CHR_ID)

ft$addRegionFeature(tbl.gwas, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()

tbl.sub.1 <- subset(tbl, atac > 0 & abs(spearmanCoeff) > 0.9 & gh)
dim(tbl.sub.1)
sort(table(tbl.sub.1$tf), decreasing=TRUE)

track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                    autoscale=TRUE)
displayTrack(igv, track)

tbl.sub.2 <- subset(tbl, atac > 0 & abs(spearmanCoeff) > 0.9 & abs(tssDistance) < 20000)
tbl.sub.2 <- subset(tbl, atac > 0 & rfScore > 0.5 & abs(tssDistance) < 20000)
dim(tbl.sub.2)


tfs <- names(sort(table(tbl.sub.2$tf), decreasing=TRUE))
for(TF in tfs){
   print(TF)
   track <- DataFrameAnnotationTrack(TF, subset(tbl.sub.2, tf==TF), color="random")
   displayTrack(igv, track)
   }
