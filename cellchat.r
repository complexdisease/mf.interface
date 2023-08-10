library(Seurat)
library(CellChat)
library(reticulate)
wdir="/path/to/rds/"
setwd(wdir)
seu=readRDS("dummy.rds")
metadata <- read.table("./metadata.tsv",sep="\t",header=T,row.names = 1)
seu@meta.data <- metadata[rownames(seu@meta.data),]
Idents(seu)="subtype"
data.input <- GetAssayData(seu, assay = "SCT", slot = "data") # normalized data matrix
labels <- Idents(seu)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
groups=c("iEVT",'eEVT','pEVT','CD14','dNK','T','cDC','DSC_A_2','DSC_B_2','aEC','vEC')
cell.use <- rownames(meta)[meta$group %in% groups] # extract the cell names from disease data

data.input <- data.input[, cell.use]
meta <- data.frame(labels = meta[cell.use,], row.names = colnames(data.input))
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP"

##LR pair
df.net <- subsetCommunication(cellchat)

groupSize=length(groups)
mat <- cellchat@net$weight[groups,groups]
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
