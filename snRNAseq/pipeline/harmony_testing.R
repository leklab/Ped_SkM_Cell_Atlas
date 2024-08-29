library(Seurat)
library(harmony)


seu_ctrl100 <- readRDS("CTRL100_doubletfiltered.rds")
seu_ctrl38 <- readRDS("CTRL38_doubletfiltered.rds")
seu_ctrl40 <- readRDS("CTRL40_doubletfiltered.rds")
seu_ctrl77 <- readRDS("CTRL77_doubletfiltered.rds")
seu_ctrl82 <- readRDS("CTRL82_doubletfiltered.rds")
seu_ctrl83 <- readRDS("CTRL83_doubletfiltered.rds")
seu_ctrl86 <- readRDS("CTRL86_doubletfiltered.rds")
seu_ctrl87 <- readRDS("CTRL87_doubletfiltered.rds")

seurat_obj <- merge(seu_ctrl100, y = c(seu_ctrl38, seu_ctrl40, seu_ctrl77, seu_ctrl82, seu_ctrl83, seu_ctrl86, seu_ctrl87), add.cell.ids = c("CTRL100", "CTRL38", "CTRL40", "CTRL77", "CTRL82", "CTRL83", "CTRL86", "CTRL87"))

table(seurat_obj$orig.ident)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)


seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", reduction.use = "pca", dims = 1:30)

# Run UMAP using the Harmony-corrected embeddings
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)

# Perform clustering on the Harmony-corrected data
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

saveRDS(seurat_obj, "harmony_combined.rds")