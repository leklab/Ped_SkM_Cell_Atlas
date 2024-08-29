#!/vast/palmer/apps/avx2/software/R/4.3.2-foss-2022b/bin/Rscript

# starsolo_seurat.R
## R script to filter and perform standard seurat workflow on STARsolo count matrices after SoupX

# Get command line arguments
args=(commandArgs(TRUE))

# Load libraries  & helper functions ----
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
# library(future) #TODO add multicore for Seurat w/ `future` for `getClusterIDs()`

suppressPackageStartupMessages(library(SoupX))

# write_sparse requirements
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(R.utils))

# source("~/DWM_utils/sc_utils/seurat_helpers/seutils.R")

SOLO_DIR = args[1] # Directory to either `Gene` or `GeneFull` matrix outputs
#NCORES = args[2]

# Helper function(s) ----
## Quick function to find clusters
getClusterIDs <- function(
    toc, 
    verbose=F
){
  seu <- CreateSeuratObject(toc)
  seu <- seu %>%
    NormalizeData(verbose=verbose)%>%
    ScaleData(verbose=verbose) %>%
    FindVariableFeatures(verbose=verbose)%>%
    RunPCA(verbose=verbose) %>%
    FindNeighbors(verbose=verbose) %>%
    FindClusters(verbose=verbose)
  
  return(
    setNames(Idents(seu), Cells(seu))
  )
}

# Read in raw and filtered matrices ----
message(paste0("Reading in SoupX matrices...\n"))
seu.data = Seurat::Read10X(paste0(SOLO_DIR,'/soupx'), gene.column=1) #droplets

seu = Seurat::CreateSeuratObject(counts = seu.data, project = "seu_test", min.cells = 3, min.features = 200)