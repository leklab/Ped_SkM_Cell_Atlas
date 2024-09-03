library(DoubletFinder)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No Seurat object provided")
}

seu_name = args[1]
#droplet_rate = args[2]

seu = readRDS(seu_name)

sweep.res.list <- paramSweep(seu, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)

#TODO
nExp_poi <- round(0.075*nrow(seu@meta.data))  ## assumed 7.5% doublet ratio- needs to be tailored to data

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


seu <- doubletFinder(seu, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seu_filtered <- subset(seu, subset = DF.classifications_0.25_0.09_255 == "Singlet")

saveRDS(seu_filtered, paste(seu_name, "_doublets.rds", sep = "")