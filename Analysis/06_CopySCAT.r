#####INITIAL ONE-TIME SETUP #####
#STEP 1: replace  "~/CopyscAT" with wherever you git cloned the repo to
#load the package
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
#parser$add_argument("--id", default=NULL,type='character',help="output_ID")
parser$add_argument("-o","--outdir",default=NULL,type='character',help="output_dir")
args <- parser$parse_args()


#use it to save references - this will create a genome chrom.sizes file, a cytobands file and a CpG file
#NOTE: updated 21-Oct-2021 with a fallback to the REST API
#generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38",tileWidth = 1e6,outputDir = "~")
DIR=args$outdir
setwd(DIR)
ID=args$id
##### REGULAR WORKFLOW #####
#initialize the environment (note: this needs to be done with every session in which you run Copy-scAT)
initialiseEnvironment(genomeFile="/mnt/data1/genome/copyscat_reference/hg38_chrom_sizes.tsv",
                      cytobandFile="/mnt/data1/genome/copyscat_reference/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="/mnt/data1/genome/copyscat_reference/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

#for this tutorial we will use the sample data included in scData
#to create your own use the process_fragment_file.py script included in the package and run it on a fragments.tsv.gz file of your choosing
#SET OUTPUT DEFAULT DIRECTORY AND NAME
setOutputFile(DIR,ID)

#PART 1: INITIAL DATA NORMALIZATION
#step 1 normalize the matrix
#USING SAMPLE DATA FROM PACKAGE
#option: if using your own file replace below with the following
scData<-readInputTable(paste0(DIR,ID,"/processed_output.tsv"))
scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)
write.table(scData_k_norm,paste0(DIR,"/",ID,"_KNorm_mat"),row.names=F,quote=F,sep="\t")

#collapse into chromosome arm level
summaryFunction<-cutAverage
scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 2,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 
#apply additional filters
scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)
#show unscaled chromosome list
graphCNVDistribution(scData_collapse,outputSuffix = ID)
#compute centers
median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)

#PART 2: ASSESSMENT OF CHROMOSOME-LEVEL CNVs 
#OPTION 1: identify chromosome-level amplifications using all cells to generate 'normal' control
#identify chromosome-level amplifications
candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=600,fakeCellSD = 0.08, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff= 0.2,medianQuantileCutoff = 0.4)
#cleanup step
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5) #= 1.5)
#final results and annotation
write.table(candidate_cnvs_clean[[1]],paste0(DIR,"/",ID,"_cell_assignment"),row.names=T,quote=F,sep="\t")
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.5,minAlteredCells=0,filterResults=FALSE,filterRange=0.8)

#PART 3: identify double minutes / amplifications
#option to compile this code
library(compiler)
dmRead<-cmpfun(identifyDoubleMinutes)
#minThreshold is a time-saving option that doesn't call changepoints on any cell with a maximum Z score less than 4 - you can adjust this to adjust sensitivity of double minute calls (note - lower value = slower)
dm_candidates<-dmRead(scData_k_norm,minCells=100,qualityCutoff2 = 100,minThreshold = 4) 
write.table(x=dm_candidates,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_dm.csv"),quote=FALSE,row.names = FALSE,sep=",")

#PART 4: assess putative LOH regions
#note: this is in beta, interpret results with caution
loh_regions<-getLOHRegions(scData_k_norm,diffThreshold = 3,lossCutoff = -0.75,minLength = 2e6,minSeg=2,targetFun=IQR,lossCutoffCells = 200,quantileLimit=0.2,cpgCutoff=100,dummyQuantile=0.6,dummyPercentile=0.4,dummySd=0.1)
write.table(x=loh_regions[[1]],file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_loss.csv"),quote=FALSE,row.names = FALSE,sep=",")

#PART 5: quantify cycling cells
#this uses the signal from a particular chromosome (we use chromosome X as typically not altered in our samples) to identify 2N versus 4N DNA content within cells
#if there is a known alteration in X in your samples, try using a different chromosome
barcodeCycling<-estimateCellCycleFraction(scData,sampName=ID,cutoff=1000)
#take max as 
write.table(barcodeCycling[order(names(barcodeCycling))]==max(barcodeCycling),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cycling_cells.tsv"),sep="\t",quote=FALSE,row.names=TRUE,col.names=FALSE)

