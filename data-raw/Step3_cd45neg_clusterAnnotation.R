## *******************************************************************
## Analyze non-immune cells: using 30PCs and annotate clusters
## Input from Immune_NonImmune_split.R
## run in 3.6.3. Seurat version 3.1.5
## *******************************************************************
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(NMF)
  library(reshape2)
})

source("source_fxns4seurat.R")

## loading Panglao database and keep major cell types
db_folder <- c("data/")
fin_db1 <- paste(db_folder,"PanglaoDB_markers_27_Mar_2020.tsv",sep="")
panglao <- read.delim(fin_db1,sep="\t",stringsAsFactors=FALSE)

listOrgan = c("Lungs","Immune system","Epithelium",
          "Connective tissue","Vasculature")
panglao = panglao[is.element(panglao$organ,listOrgan),]

listSelCT = c("Pulmonary alveolar type I cells",
              "Pulmonary alveolar type II cells",
              "Clara cells","Airway epithelial cells",
              "Airway goblet cells","Ciliated cells",
              "Ionocytes","Fibroblasts",
              "Endothelial cells")

tmpsel  = grepl("Hs",panglao$species) & 
          is.element(panglao$cell.type,listSelCT)
db_sub = panglao[tmpsel,2:3]
colnames(db_sub)=c("Symbol","cellType")

## Step 1: load Seurat object
fin1 <- c("GGO_seuInte_cd45neg.rds")
cd45neg <-readRDS(fin1)

## Step 2: scale and dimentional reduction
ndim=30  #*****
cd45neg <- FindNeighbors(cd45neg, dims = 1:ndim)
cd45neg <- FindClusters(cd45neg, resolution = 1)
cd45neg <- RunUMAP(cd45neg,dims = 1:ndim)

# input for cluster annotation
target_seu <- cd45neg   

## Step 3: cluster annotation
## Step 3.1: enrichment approach
## Step 3.1.1: find cluster markers
allClust_mks <- FindAllMarkers(target_seu, only.pos = TRUE, 
         min.pct = 0.25, logfc.threshold = 0.25)

## export data
##write.table(allClust_mks, c("GGO_cd45neg_markers.txt"),
##          quote = F, row.names = F)

## Step 3.1.2: enrichment approach by using panglao database
enrichedData <- enrichScoreCalc_mt(target_seu,
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

## Step 3.2: using cannonical markers from either literature (FACS markers) 
## or overlap genes based on enrichment analysis
gList = c("COL1A1","COL1A2", #fibro
         "PECAM1","ECSCR",   #endo 
          "EPCAM","NKX2-1","SFTPC","SFTPA1",  ##AT2 or Clara
          "SCGB1A1","SCGB3A1",  #Clara or BASC
          "MUC5AC","MUC16",   #Goblet
          "PDPN",   #AT1 
          "TFCP2L1",  #ionocyte, 
          "FOXJ1")    #ciliated

# Get (arithmetic) average of markers: non-log scale
tmp_averages <- AverageExpression(target_seu,
        assays="RNA",
        features=gList)$RNA
tmp_averages <- t(tmp_averages)
tmp_zscores <- scale(log(tmp_averages+1))   ##high dynamic rank --> log transform to norm distribution
mks_clust_avg <- cbind(rownames(tmp_zscores),tmp_zscores)

## Step 3.3: export output from enrichment and marker expression
Ftab_out <- cbind(mks_clust_avg,Ftab_candidate)
fout <- c("data/GGO_cd45neg_pred.txt")
write.table(Ftab_out, fout, sep="\t",
          quote = F, row.names = F)

## Plot z-scores of marker expression --> (+) threshold
tmp4plot <- melt(tmp_zscores)
names(tmp4plot) <- c("clustID","Gene","Expr")
tmpOrd <- order(tmp4plot$Gene,tmp4plot$Expr,decreasing = FALSE)
tmp4plot <- tmp4plot[tmpOrd,]
tmp4plot$index <- rep(c(1:nrow(tmp_zscores)),times=ncol(tmp_zscores))

p1 <- ggplot(tmp4plot, aes(x=index,y=Expr)) + 
    geom_point()+
    facet_wrap(~Gene, scale="free",ncol=5)+
    theme_bw()

## Heatmap of enrichment scores
p2 <- aheatmap(enrichedData$score_Mt,scale="col")

## Step 3.4: Manually inspect the 3.3 output and assign cell type:
## (1) if difference between 1st and 2nd highest score > 1, 
##         assign cell type associated with 1st highest score
## (2) else
##         assigne based on cell type expressing canonical markers
##         AT2 = EPCAM and NKX2-1 
##         BASC = SCGB1A1 and SCGB1A1 and PECAM (intermediated level)
##         Clara = EPCAM, NKX2-1, SCGB1A1,SCGB1A1
## Create matrix of clusters X cell_Type for next step

## Step 4: Update cell types after step 3.4
fin3 <- c("example/GGO_cd45neg_CID_annot.txt")
cellAnnotInfo <- read.delim(fin3,sep="\t",header=T,
        stringsAsFactors=FALSE)

cellID_vec <- cellAnnotInfo[,2]
names(cellID_vec) <- cellAnnotInfo[,1]

nonIm_level1 <- c("AT1","AT2","AT2_loSFTPC","Clara","Ciliated",
                  "endo","fibro","BASC")
nonIm_level2 <- c("AT1","AT2","Clara","Ciliated",
                  "endo","fibro", "BASC")

# add cell type annotation
tmpClust = Idents(cd45neg)
cd45neg <- AddMetaData(cd45neg,
           metadata=tmpClust,col.name='cd45neg_clustID')

cd45neg <- RenameIdents(cd45neg,cellID_vec)
tmp = Idents(cd45neg)
tmp2 = factor(tmp,levels=nonIm_level1,
           labels=nonIm_level1)
cd45neg <- AddMetaData(cd45neg,tmp2,
                   col.name = 'cd45neg_cellID')

# add histology as binary levels
tmpBin = ifelse(cd45neg@meta.data$histology=="Normal",
            "Normal","GGO_Tu")
tmpBin = factor(tmpBin,levels=c("Normal","GGO_Tu"),
           labels=c("Normal","GGO_Tu"))
cd45neg <- AddMetaData(cd45neg,tmpBin,
                   col.name = 'binaryHis')

saveRDS(cd45neg,file=c("GGO_cd45neg_seu.rds"))

## Step 5 (Option): export data of a specific cell type
## Below is an example of exporting AT2 cells for 
## analyses (e.g. DEG analysis)
Idents(cd45neg) <- 'cd45neg_cellID'
selID <- c("AT2","AT2_loSFTPC")   #*****
target_seu <- subset(cd45neg, idents = selID, invert = FALSE) #*****
tmp <- Idents(target_seu)
target_seu <- AddMetaData(target_seu,
           metadata=tmp,col.name='cd45neg_oneAT2')
Idents(target_seu) <- 'cd45neg_clustID'    ##convert back to cluster ID
saveRDS(target_seu,file=c("GGO_AT2_seu.rds"))

## a subset of fibroblast (<100MB size)
Idents(cd45neg) <- 'cd45neg_cellID'
selID <- c("fibro")   #*****
target_seu <- subset(cd45neg, idents = selID, invert = FALSE) #*****
Idents(target_seu) <- 'cd45neg_clustID'    ##convert back to cluster ID
sel_cID <- c(9,11)
target_seu <- subset(target_seu, idents = sel_cID, invert = FALSE) #*****
saveRDS(target_seu,file=c("GGO_fibro_seu_gitTest.rds"))
