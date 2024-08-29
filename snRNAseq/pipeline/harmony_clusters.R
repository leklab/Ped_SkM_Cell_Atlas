library(Seurat)
library(harmony)


harmony_combined <- readRDS("harmony_nosoup.rds")

harmony_combined[["RNA"]] <- JoinLayers(harmony_combined[["RNA"]])

Idents(harmony_combined) <- "seurat_clusters"


# Find markers for all clusters
all_markers <- FindAllMarkers(harmony_combined, only.pos = TRUE)

# View the top markers for each cluster
head(all_markers)

write.table(all_markers, file = "all_markers_harmony_nosoup.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
