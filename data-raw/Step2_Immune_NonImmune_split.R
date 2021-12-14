## **************************************************************
## Splitting immune and non-immune cells by canonical markers
## Input from DataPreProcessing_BatchAdjustment.R
## Require reference gene sets (e.g. Panglao database)
## run in 3.6.3. Seurat version 3.1.5
## **************************************************************
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(NMF)
})

source("R/source_fxns4seurat.R")

## loading Panglao database and keep major cell types
db_folder <- c("data/")
fin_db1 <- paste(db_folder,"PanglaoDB_markers_27_Mar_2020.tsv",sep="")
panglao <- read.delim(fin_db1,sep="\t",stringsAsFactors=FALSE)

listOrgan = c("Lungs","Immune system","Epithelium",
          "Connective tissue","Vasculature")
panglao = panglao[is.element(panglao$organ,listOrgan),]


listSelCT = c("T cells","B cells","Dendritic cells",
                "Macrophages","Alveolar macrophages",
                "Mast cells","Monocytes","Neutrophils",
                "NK cells","Plasma cells",
                "Plasmacytoid dendritic cells",
                "Pulmonary alveolar type I cells",
                "Pulmonary alveolar type II cells",
                "Clara cells","Airway epithelial cells",
                 "Airway goblet cells","Ciliated cells",
                 "Ionocytes","Fibroblasts","Stromal cells",
                 "Endothelial cells")

tmpsel  = grepl("Hs",panglao$species) & 
          is.element(panglao$cell.type,listSelCT)
db_sub = panglao[tmpsel,2:3]
colnames(db_sub)=c("Symbol","cellType")

## Step 1: load Seurat object
fin1 <- c("GGO_BatchAdj.rds")
allSam.integrated <-readRDS(fin1)

DefaultAssay(allSam.integrated) <- "integrated"

## Step 2: scale and dimentional reduction
allSam.integrated <- ScaleData(allSam.integrated, verbose = FALSE)
allSam.integrated <- RunPCA(allSam.integrated, npcs = 100, verbose = FALSE)
allSam.integrated <- FindNeighbors(allSam.integrated, dims = 1:50)
allSam.integrated <- FindClusters(allSam.integrated, resolution = 1)
allSam.integrated <- RunUMAP(allSam.integrated, dims = 1:49)
           
## Step 3: cluster annotation
## Step 3.1: enrichment approach
## Step 3.1.1: find cluster markers
allClust_mks <- FindAllMarkers(allSam.integrated, only.pos = TRUE, 
         min.pct = 0.25, logfc.threshold = 0.25)

## export data for record
##fin_4mks <- c("GGO_allCells_markers.txt")
##tmpOut <- cbind(rownames(tmp_tab),tmp_tab)
##write.table(tmpOut,fin_4mks,sep="\t",quote=F,row.names=F)


## Step 3.1.2: enrichment approach by using panglao database
enrichedData <- enrichScoreCalc_mt(allSam.integrated,
                    allClust_mks,db_sub)

## Find the top two highest scores and associated cell types for each cluster
nclusters <- nrow(enrichedData$score_Mt)
Ftab_candidate <- matrix(NA,nclusters,4)

for (i in c(1:nclusters)){
   tmpi_hi <- order(enrichedData$score_Mt[i,],decreasing = TRUE)[c(1,2)]
   tmpscore_hi <- enrichedData$score_Mt[i,tmpi_hi]
   tmpCell_hi <- colnames(enrichedData$score_Mt)[tmpi_hi]
   Ftab_candidate[i,] <- c(tmpscore_hi,tmpCell_hi)
}

colnames(Ftab_candidate) <- c("score1","score2","cellType1","cellType2")
rownames(Ftab_candidate) <- rownames(enrichedData$score_Mt)

## Step 3.2: using cannonical markers from literature
gList = c("PTPRC",    #CD45 for immune
          "CD3D",     #Tcells
          "CD79A",    #B and plasma
          "CD68",     #myeloid
          "COL1A1",   #fibroblast
          "PECAM1",   #CD31 for endothelial
          "EPCAM","SCGB1A1")  #lungg epithelial marker    

# Get (arithmetic) average of markers: non-log scale
tmp_averages <- AverageExpression(allSam.integrated,
        assays="RNA",
        features=gList)$RNA
tmp_averages <- t(tmp_averages)
mks_clust_avg <- cbind(rownames(tmp_averages),tmp_averages)

## Step 3.3: export output from enrichment and marker expression
##Ftab_out <- cbind(mks_clust_avg,Ftab_candidate)
##fout <- c("GGO_allClusters_pred_040921.txt")
##write.table(Ftab_out, fout, sep="\t",
##          quote = F, row.names = F)

## Step 3.4: Manually inspect the 3.3 output and categorize into:
## endothelial --> fibroblast --> epithelial, including secteroty cells 
## --> immune based on: 
## (1) if difference between 1st and 2nd highest score > 1, 
##         assign cell type associated with 1st highest score
## (2) else
##         assigne based on cell type whose marker expression is highest
## Create matrix of clusters X cell_Type for splitting immune vs non-immune

## Step 4: Splitting data into immune and non-immune cells
## load anotation matrix created after Step 3.3
finAnnot <- c("example/allCells_cID_annot.txt")
tmp <- read.delim(finAnnot,sep="\t",head=T,
            stringsAsFactors=FALSE)

cellID_vec <- tmp[,2]
names(cellID_vec) <- tmp[,1]

## check available cell attributions
colnames(allSam.integrated@meta.data)

## add cell type annotation
tmpClust = Idents(allSam.integrated)
allSam.integrated <- AddMetaData(allSam.integrated,
           metadata=tmpClust,col.name='clustID')

allSam.integrated <- RenameIdents(allSam.integrated,cellID_vec)
tmp = Idents(allSam.integrated)
allSam.integrated <- AddMetaData(allSam.integrated,tmp,
                   col.name = 'cellIDAll')

# Split data into 2 groups: immune vs. non-immune
saveRDS(allSam.integratedt,file=c("GGO_seuInte_all.rds"))

cd45pos_seurat <- subset(allSam.integrated,idents='immune')
saveRDS(cd45pos_seurat,file=c("GGO_seuInte_cd45pos.rds"))

cd45neg_seurat <- subset(allSam.integrated,idents='immune',invert=TRUE)
saveRDS(cd45neg_seurat,file=c("GGO_seuInte_cd45neg.rds"))


