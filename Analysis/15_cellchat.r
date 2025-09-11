library(Seurat)
library(CellChat)
library(reticulate)
## loading output data from stereo-seq
df=read.table("./spatial_celltype_anno.tsv",sep="\t",row.names = 1,head =T)
data.input= t(readMM('/path/to/stereoseq_gex.mtx'))
meta=read.table("stereoseq_obs.tsv",sep="\t",row.names = 1,head=T)
vars=read.table("stereoseq_var.tsv",sep="\t",row.names = 1,head=T)
rownames(data.input)=rownames(vars)
colnames(data.input)=rownames(meta)

wdir="/path/to/rds/"
setwd(wdir)
spatial.factors = data.frame(ratio = 10/20, tol = 5)

meta$labels <- meta[["subclass"]] %>% as.character()
meta[meta$subclass=='EVT','labels']=meta[meta$subclass=='EVT','subtype'] %>% as.character()
spatial.locs=meta[,c('shift_x','shift_y')]
color.use <- scPalette(length(table(meta$labels))); names(color.use) <- names(table(meta$labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = list(spot.diameter=10,spot=1.8))

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 8000000 * 1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean",
                                distance.use = TRUE, interaction.length = 250, scale.distance = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

target=c("immune", "DSC", "Epi", "FB", "fVEC", "HB", "mVEC","PV","SCT","VCT",'eEVT','iEVT')
groupSize=length(target)
mat <- cellchat@net$weight
mat=mat[target,target]
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  #mat=mat[target,target]
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  #mat2[i,i]=0
  #scale_up=sort(as.vector(mat[,i]), decreasing = TRUE)[2]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.8,edge.weight.max = max(mat2), title.name = rownames(mat)[i])
}
##permutation
pval_mat= matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
testA='iEVT'

target=c("dNK/T/dM", "DSC", "Epi", "FB", "fVEC", "HB", "mVEC","PV","SCT","VCT",'eEVT','iEVT')
for (testB in target){
observed_weight <- cellchat@net$weight[testA, testB]

# Manually permute cell type labels and recompute weights
n_perm <- 1000
perm_weights <- list()

for (i in 1:n_perm) {
  permuted <- cellchat
  permuted@idents <- sample(cellchat@idents)  # shuffle cell type labels
  permuted <- computeCommunProb(permuted)
  permuted <- computeCommunProbPathway(permuted)
  permuted <- aggregateNet(permuted)
  perm_weights[[i]] <- permuted@net$weight
}
p_value <- mean(perm_weights >= observed_weight)
pval_mat[testA,testB]=p_value
}
observed_weight=cellchat@net$weight

tmp=array(0, dim = c(nrow(perm_weights[[1]]), ncol(perm_weights[[1]]), length(perm_weights)))
for (i in 1:100){
  tmp[,,i]=perm_mat[,,i]>=observed_weight
}
pval_mat=matrix(0,ncol = ncol(perm_weights[[1]]),nrow = nrow(perm_weights[[1]]))
for (i in 1:100){
  pval_mat=pval_mat+tmp[,,i]
}



                                          
