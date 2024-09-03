library(Seurat)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No SoupX directory provided")
}

SOLO_DIR = args[1]
sample_name = args[2]

# Read in SoupX matrices ----
message(paste0("Reading in SoupX matrices...\n"))
seu.data = Seurat::Read10X(paste0(SOLO_DIR,'/soupx'), gene.column=1) 

seu = Seurat::CreateSeuratObject(counts = seu.data, project = sample_name, min.cells = 3, min.features = 200)

#Standard Seurat QC/Filtering Metrics

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 5)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

saveRDS(seu, file = paste(sample_name, ".rds", sep = ""))

#Finding markers
#Idents(seu) <- "seurat_clusters"


# Find markers for all clusters
#all_markers <- FindAllMarkers(seu, only.pos = TRUE)

# View the top markers for each cluster
#head(all_markers)

#write.table(all_markers, file = paste(sample_name, "all_markers.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
