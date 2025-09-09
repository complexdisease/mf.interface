## Rscript of converting gene expression seurat object to h5ad 
library(Seurat)
library(SeuratData)
library(SeuratDisk)
seu=readRDS("/path/to/the/integrated_gex.rds")
SaveH5Seurat(seu, filename = "/path/to/the/integrated_gex.h5Seurat")
Convert("/path/to/the/integrated_gex.h5Seurat", dest = "h5ad")
