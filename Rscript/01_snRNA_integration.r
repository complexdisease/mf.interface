library(Seurat)
library(future)
library(dplyr)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 200000 * 1024^2) # for 50 Gb RAM

dir="/path/to/cellranger/output/"
setwd(wdir)
files=list.files(dir)
bc=list()
metric_list=list()
mat_list=list()
assay_list=list()
obj_list=list()

for (i in 1:length(files)){
        bc[[i]]=read.table(gzfile(paste0(dir,files[i],"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")),header=F)
        metric_list[[i]]=read.csv(paste0(dir,files[i],"/outs/per_barcode_metrics.csv"),header=T,row.names=1)
	metric_list[[i]]=metric_list[[i]][bc[[i]]$V1,]
}

for (i in 1:length(files)){counts <- Read10X_h5(paste0(dir,files[i],"/outs/filtered_feature_bc_matrix.h5"))
	obj_list[[i]] <- CreateSeuratObject(counts = counts$`Gene Expression`, project = files[i], min.cells = 3, min.features = 10,assay="RNA",meta.data=metric_list[[i]])
	obj_list[[i]]$dataset=files[i]
	obj_list[[i]][["percent.mt"]] <- PercentageFeatureSet(obj_list[[i]], pattern = "^MT-")
	obj_list[[i]]@meta.data$barcode=rownames(obj_list[[i]]@meta.data)
	obj_list[[i]]=RenameCells(obj_list[[i]],add.cell.id=files[i])
}

obj_list <- lapply(X = obj_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20,reference = 2)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30,k.weight=20)
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = 0.5)
saveRDS(combined.sct,"integrated_gex.rds")

