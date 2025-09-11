library(biomaRt)
library(dplyr)
library(Signac)
library(patchwork)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RColorBrewer)
setwd("/path/to/workdir/")
seu=readRDS("/path/to/integrated_atac_filtered.rds")
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("external_gene_name",'entrezgene_id'), values=names(genes),filters ='entrezgene_id', mart = mart)
names(genes) <- bm$external_gene_name[match(genes$gene_id,bm$entrezgene_id)]

colors=colorRampPalette(brewer.pal(8, "Spectral"))(15)
options(ggplot2.discrete.color=colors)
PR <- promoters(genes,downstream=400,upstream = 2000)
a@meta.data$subclass=factor(a@meta.data$subclass,levels = c('VCT',"EVT","SCT","STM","FB","PV","Epi",
                                                "Endo","B","T","dNK","M","HB","cDC","Erythrocyte"))
g1=c("TENM3",'HLA-G',"CYP19A1",'DKK1','BMP5',"MEOX2",'MECOM','CDH5','BCL11A','CD3D','GNLY','MS4A7',"CD74",'HBB')
Idents(a)="subclass"
p=list()
for (i in 1:length(g1)){
p[[i]]=CoveragePlot(object = a,region = PR[g1[i]],sep = c(":", "-"),group.by="subclass",extend.upstream = 0,
extend.downstream = 500,ncol = 1,peaks = F) + scale_fill_manual(values=colors)}
p1 <- cowplot::plot_grid(plotlist=p,ncol = 14)
