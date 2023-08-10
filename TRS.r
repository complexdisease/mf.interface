library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc = "/tank/data2/cw/miniconda/miniconda3/envs/r4/lib/R/library")
library(igraph)
seu=readRDS("/path/to/rds/")
setwd("/working/directories/")
set.seed(9527)
peaks=brain[['ATAC']]@ranges
peaks=keepSeqlevels(peaks,seqlevels(peaks)[1:22],pruning.model='coarse')
SE <- SummarizedExperiment(assays = list(counts = seu[['ATAC']]@counts[1:length(peaks),]),
                           rowData = peaks, 
                           colData = seu@meta.data)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38)
SE_bg <- getBackgroundPeaks(SE, niterations=200)
ukbb <- importBedScore(rowRanges(SE), "sPTB.PP.abf.hg38.bed", colidx = 5)
SE_DEV <- computeWeightedDeviations(object = SE, weights = ukbb, background_peaks = SE_bg)
z_score_mat <- data.frame(SummarizedExperiment::colData(SE), 
                          z_score=t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c())
seed_idx <- seedindex(z_score = z_score_mat$z_score,
                      percent_cut = 0.05)
scale_factor <- cal_scalefactor(z_score = z_score_mat$z_score, 
                                percent_cut = 0.01)
peak_by_cell_mat <- SummarizedExperiment::assay(SE)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
lsi_mat=Embeddings(object = seu[["integrated_lsi"]])[,1:30]
mutualknn30 <- getmutualknn(lsimat = lsi_mat2, 
                            num_k = 30)
np_score <- randomWalk_sparse(intM = mutualknn30, 
                              queryCells = rownames(mutualknn30)[seed_idx], 
                              gamma = 0.05)
omit_idx <- np_score==0
sum(omit_idx)
mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- capOutlierQuantile(x = np_score, 
                          q_ceiling = 0.95) |> max_min_scale()
TRS <- TRS * scale_factor
mono_mat <- data.frame(z_score_mat[!omit_idx, ], 
                       seed_idx[!omit_idx], 
                       np_score, 
                       TRS)
write.table(mono_mat,"sPTB.TRS",sep="\t",row.names = T)


mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, 
                                 seed_idx=mono_mat$seed_idx, 
                                 topseed_npscore=mono_mat$np_score, 
                                 permutation_times=1000, # Increase to >=1000 in practice
                                 true_cell_significance=0.05, 
                                 rda_output=FALSE, 
                                 # mycores=8,# Increase in practice
                                 rw_gamma=0.05)
mono_mat2 <- data.frame(mono_mat, mono_permu)

##enriched cells
mono_mat2 |>
  dplyr::group_by(color) |> 
  dplyr::summarise(enriched_cell=sum(true_cell_top_idx)) |> 
  ggplot(aes(x=color, y=enriched_cell, fill=color)) + 
  geom_bar(stat="identity") + 
  theme_classic() 
mono_mat2$rev_true_cell_top_idx <- !mono_mat2$true_cell_top_idx
## depleted cells
mono_mat2 |>
  dplyr::group_by(color) |> 
  dplyr::summarise(depleted_cell=sum(rev_true_cell_top_idx)) |> 
  ggplot(aes(x=color, y=depleted_cell, fill=color)) + 
  geom_bar(stat="identity") + 
  theme_classic()
