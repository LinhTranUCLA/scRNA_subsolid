## Preprocessing data: Combine samples, filter and batch adjustment
## run in 3.6.3. Seurat version 3.1.5
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

## Step 0: specify working folder and data location 
workFolder = c("scRNA_nodules/Ranalysis_Gencode34/")
setwd(workFolder)

inputFolder= c("scRNA_nodules/CountMatrices/")


## Step 1: load data
## Step 1a: specify input files
sampID_vec = c("Case1_Normal","Case1_GGO",
               "Case2_Normal","Case2_GGO","Case2_Tumor",
               "Case3_Normal","Case3_GGO",
               "Case4_Normal","Case4_GGO1","Case4_GGO2","Case4_Tumor",
               "Case5_Normal","Case5_GGO",
               "Case6_Normal","Case6_GGO") ##***

sampID_short = c("C1N","C1G",
                 "C2N","C2G","C2T",
                 "C3N","C3G",
                 "C4N","C4G1","C4G2","C4T",
                 "C5N","C5G",
                 "C6N","C6G") ##***

folderID_vec = sapply(sampID_vec,function(x)
     paste0(unlist(strsplit(x,split="_")),collapse="",sep=""))
filterFolder = c("/filtered_feature_bc_matrix/")
refOption = c("_singlePAR/outs")

allPaths <- paste(inputFolder,folderID_vec,refOption,filterFolder,sep="")

# Step 1b: import 10X data and create seurat objects
nsamples = length(allPaths)
allDataRaw = list()
for (sample.Idx in c(1:nsamples)){
   tmp1 <- Read10X(data.dir = allPaths[sample.Idx])
   tmp1_seu <- CreateSeuratObject(counts = tmp1,
        project = sampID_vec[sample.Idx])
   allDataRaw[[sample.Idx]] <- tmp1_seu
  
}

## Step 1c: Combine seurat object based on batches
batch1 <- merge(x=allDataRaw[[1]],y=allDataRaw[[2]],
              add.cell.ids=sampID_short[1:2],
              project="Pat1")
batch2 <- merge(x=allDataRaw[[3]],y=c(allDataRaw[[4]],allDataRaw[[5]]),
              add.cell.ids=sampID_short[3:5],
              project="Pat2")
batch3 <- merge(x=allDataRaw[[6]],y=c(allDataRaw[[7]],allDataRaw[[8]],
              allDataRaw[[9]],allDataRaw[[10]],allDataRaw[[11]],
              allDataRaw[[12]],allDataRaw[[13]],allDataRaw[[14]],
              allDataRaw[[15]]),add.cell.ids=sampID_short[6:15],
              project="Pat3_6")

## Step 2: QC and filter
batch1[["percent.mt"]] <- PercentageFeatureSet(batch1, pattern = "^MT-")
batch2[["percent.mt"]] <- PercentageFeatureSet(batch2, pattern = "^MT-")
batch3[["percent.mt"]] <- PercentageFeatureSet(batch3, pattern = "^MT-")

## Filtered out "bad" cells
batch1_2 <- subset(batch1, subset = nFeature_RNA > 475 & percent.mt < 17)
batch2_2 <- subset(batch2, subset = nFeature_RNA > 475 & percent.mt < 17)
batch3_2 <- subset(batch3, subset = nFeature_RNA > 475 & percent.mt < 17)

## Step 3: Normalize and Find HVGs
allSam.list <- list(batch1_2,batch2_2,batch3_2)

for (i in 1:length(allSam.list)) {
    allSam.list[[i]] <- NormalizeData(allSam.list[[i]], verbose = FALSE)
    allSam.list[[i]] <- FindVariableFeatures(allSam.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

## Step 4: Integrate data by seurat pipeline
allSam.anchors <- FindIntegrationAnchors(object.list = allSam.list, 
                     dims = 1:100)
allSam.integrated <- IntegrateData(anchorset = allSam.anchors, dims = 1:100)

fout1 <- c("GGO_BatchAdj.rds")
saveRDS(allSam.integrated,file=fout1)




