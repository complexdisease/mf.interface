library(GenomicRanges)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

## set.seed(1234)
## function of converting bed files to Granges 
bed_to_granges <- function(file){
df <- read.table(file,
header=F,
stringsAsFactors=F)
if(length(df) > 6){
df <- df[,-c(7:length(df))]
}
if(length(df)<3){
stop("File has less than 3 columns")
}
header <- c('chr','start','end','id','score','strand')
names(df) <- header[1:length(names(df))]
if('strand' %in% colnames(df)){
df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
}
if(length(df)==3){
gr <- with(df, GRanges(chr, IRanges(start, end)))
} else if (length(df)==4){
gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
} else if (length(df)==5){
gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
} else if (length(df)==6){
gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
}
return(gr)
}

atac=readRDS('/path/to/integrated_ATAC_filtered.rds')  ## Load post-QC snATAC-seq obj
enhancers=bed_to_granges("/path/to/Fantom5.hg38.enhancers.bed"). ## Load the downloaded FANTOM5 enhancers
enhancers <- keepStandardChromosomes(enhancers, pruning.mode = "coarse")
enhancers<- subsetByOverlaps(x = enhancers, ranges = blacklist_hg38_unified, invert = TRUE)
counts <- FeatureMatrix(fragments = Fragments(atac),features = enhancers,cells = colnames(atac))

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) ## Genome annotation
seqlevelsStyle(annotation) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(enhancer_obj) <- annotations

chrom_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), genome = 'hg38',fragments = Fragments(atac)
min.cells = 0, min.features=0)
enhancer_obj <- CreateSeuratObject(counts = chrom_assay,assay = "enhancers")
enhancer_obj@meta.data=atac@meta.data[row.names(enhancer_obj)]
SaveH5Seurat(enhancer_obj, filename = "F5_enhancer.h5Seurat")
Convert("F5_enhancer.h5Seurat", dest = "h5ad")        ## Convert to h5ad format
