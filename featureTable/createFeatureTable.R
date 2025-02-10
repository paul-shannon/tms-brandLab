library(trena)
targetGene <- "KLF1"

ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg38")
tbl.fimo <- get(load("../fimo/KLF1-20k-em5/fimo-5.RData"))
dim(tbl.fimo)  # 40664     9
ft$setFundamentalRegions(tbl.fimo)
mtx <- get(load("brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"))
dim(mtx) # 27171 28
fivenum(mtx["KLF1",])   # 7.89: 11.76

candidateTFs <- sort(unique(intersect(tbl.fimo$tf, rownames(mtx))))
length(candidateTFs) # 586
deleters <- grep("^KLF1$", candidateTFs)
if(length(deleters) > 0)
    candidateTFs <- candidateTFs[-deleters]
length(candidateTFs) # 585

solver <- EnsembleSolver(mtx,
                         targetGene=targetGene,
                         candidateRegulators=candidateTFs,
                         solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                         geneCutoff=1.0)
tbl.out <- suppressWarnings(run(solver))
new.order <- order(abs(tbl.out$spearmanCoeff), decreasing=TRUE)
tbl.out <- tbl.out[new.order,]
dim(tbl.out)  # 577 7

ft$addGeneFeature(tbl.out[, c("gene", "rfScore")],
                  feature.name="rfScore", default.value=0)
ft$addGeneFeature(tbl.out[, c("gene", "spearmanCoeff")],
                  feature.name="spearmanCoeff", default.value=0)
tbl <- ft$getTable()

tbl.atac <- get(load("tbl.atacMerged.RData"))
dim(tbl.atac)
head(tbl.atac)

feature.guide <- list(atac="width")
default.values <- list(atac=0)
ft$addRegionFeature(tbl.atac, feature.genome="hg38", feature.guide, default.values)
tbl <- ft$getTable()

tbl.gwas.raw <- get(load("tbl.gwas-vuckovic.RData"))
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
