## run in 3.6.3. Seurat version 3.1.5
## Preprocessing data: Combine filtered count matrices, filter and batch adjustment
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(RColorBrewer)
  library(ggplot2)
  library(NMF)
  library(gplots)
  library(lmerTest)
  library(limma)
  library(MAST)
  library(edgeR)
  library(reshape)
  library(iTALK)
  library(ggrepel)
})

source("R/source_fxns4seurat.R")

## main script
## Step 1: Import cell ranger output
batch1 <- importCellRangerCntMts(c("CellRangeInputDirs_batch1.txt"))
batch2 <- importCellRangerCntMts(c("CellRangeInputDirs_batch2.txt"))
batch3 <- importCellRangerCntMts(c("CellRangeInputDirs_batch3.txt"))

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




