library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

setwd('/path/to/workdir/')
seu <- readRDS("/path/to/rds")
DefaultAssay(seu) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(seu)) %in% main.chroms)
seu <- seu[keep.peaks, ]
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
seu <- AddMotifs(
  object = seu,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

##Motif activity by chromVar
seu <- RunChromVAR(
  object = seu,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
# TF Footprinting by Signac
features=c("TP63","GATA3","SPI1", "RUNX2")
seu <- Footprint(
  object = seu,
  motif.name = features,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

p2 <- PlotFootprint(seu, group.by='subclass',idents = c("VCT","EVT","SCT","M","dNK","Epi"), features = features)
