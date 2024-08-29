# BCH_analysis_pipeline
Preliminary workflow analysis for BCH single nuclei data

File names here are hard coded in, this is just to serve as a temporary reference to what is being run in the pipeline. I will clean this up and make it easier to run in the future.

Workflow 
1. Run starsolo_soupx.R on GeneFull/ or Gene/ directory from STARsolo otuput to remove ambient mRNAs
2. Run through seurat_testing.Rmd with output soupx/ directory from step 1 (standard seurat workflow)
3. Run DoubletFinder on seurat object generated from step 2 to identify and filter doublets
4. Run harmony_testing.R to produce a seurat object with all samples integrated and recluster based on new seurat object
