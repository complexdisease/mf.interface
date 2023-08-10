library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 200000 * 1024^2) # for 200 Gb RAM
dir="/path/to/cellranger/output/"
wdir="/path/to/workdir/"
setwd(wdir)
files=list.files(dir)
bc=list()
metric_list=list()
fragment_list=list()
mat_list=list()
assay_list=list()
obj_list=list()
peaks=list()
gr=list()

for ( i in 1:length(files)){
	peaks[[i]]=read.table(file=paste0(dir,files[i],"/outs/atac_peaks.bed"),col.names=c("chr", "start", "end"))
	gr[[i]]=makeGRangesFromDataFrame(peaks[[i]])
}

combined.peaks=reduce(do.call(c,gr))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

for (i in 1:length(files)){
        bc[[i]]=read.table(gzfile(paste0(dir,files[i],"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")),header=F)

        metric_list[[i]]=read.csv(paste0(dir,files[i],"/outs/per_barcode_metrics.csv"),header=T,row.names=1)
        metric_list[[i]]=metric_list[[i]][bc[[i]]$V1,]
}

for (i in 1:length(files)){
        fragment_list[[i]]=CreateFragmentObject(path=paste0(dir,files[i],"/outs/atac_fragments.tsv.gz"),cells = bc[[i]]$V1)
        mat_list[[i]]=FeatureMatrix(fragments = fragment_list[[i]],features = combined.peaks,cells = bc[[i]]$V1)
        assay_list[[i]] <- CreateChromatinAssay(mat_list[[i]], fragments = fragment_list[[i]])
        obj_list[[i]]=CreateSeuratObject(assay_list[[i]], assay = "ATAC",meta.data=metric_list[[i]])
}

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

for (i in 1: length(files)){
	Annotation(obj_list[[i]]) <- annotations
	obj_list[[i]] <- FindTopFeatures(obj_list[[i]], min.cutoff = 10)
        obj_list[[i]] <- RunTFIDF(obj_list[[i]])
        obj_list[[i]] <- RunSVD(obj_list[[i]])
	obj_list[[i]]<- NucleosomeSignal(object = obj_list[[i]])
	obj_list[[i]]<- TSSEnrichment(object = obj_list[[i]],fast=F)
	obj_list[[i]]$pct_reads_in_peaks <- obj_list[[i]]$atac_peak_region_fragments / obj_list[[i]]$atac_fragments * 100
	#obj_list[[i]]$blacklist_ratio <- obj_list[[i]]$blacklist_region_fragments / obj_list[[i]]$peak_region_fragments
        obj_list[[i]]@meta.data$dataset=files[i]
	obj_list[[i]]@meta.data$BC=rownames(obj_list[[i]]@meta.data)
	obj_list[[i]]=RenameCells(obj_list[[i]],new.names=paste(obj_list[[i]]@meta.data$dataset,obj_list[[i]]@meta.data$BC,sep="_"))
}
combined <- merge(x = obj_list[[1]],y = obj_list[2:length(files)])
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunTFIDF(combined)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:30)
integration.anchors <- FindIntegrationAnchors(object.list = obj_list,anchor.features = rownames(combined),reduction = "rlsi",dims = 2:30)
integrated <- IntegrateEmbeddings(anchorset = integration.anchors,reductions = combined[["lsi"]],new.reduction.name = "integrated_lsi",dims.to.integrate = 1:30,k.weights=100)
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:50)
saveRDS(integrated,"integrated_ATAC.rds")

