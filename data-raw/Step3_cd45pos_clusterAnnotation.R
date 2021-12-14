## *******************************************************************
## Analyze immune cells: using 30PCs and annotate clusters
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

listSelCT = c("T cells","B cells","Dendritic cells",
              "Macrophages","Alveolar macrophages",
              "Mast cells","Monocytes","Neutrophils",
              "NK cells","Plasma cells",
              "Plasmacytoid dendritic cells",
              "Gamma delta T cells")

tmpsel  = grepl("Hs",panglao$species) & 
          is.element(panglao$cell.type,listSelCT)
db_sub = panglao[tmpsel,2:3]
colnames(db_sub)=c("Symbol","cellType")

## Step 1: load Seurat object
fin1 <- c("GGO_seuInte_cd45pos.rds")
cd45pos <-readRDS(fin1)


## Step 2: scale and dimentional reduction
ndim=30  #*****
cd45pos <- FindNeighbors(cd45pos, dims = 1:ndim)
cd45pos <- FindClusters(cd45pos, resolution = 1)
cd45pos <- RunUMAP(cd45pos,dims = 1:ndim)

# input for cluster annotation
target_seu <- cd45pos  

## Step 3: cluster annotation
## Step 3.1: enrichment approach
## Step 3.1.1: find cluster markers
allClust_mks <- FindAllMarkers(target_seu, only.pos = TRUE, 
         min.pct = 0.25, logfc.threshold = 0.25)

## export data
##write.table(allClust_mks, c("GGO_cd45pos_markers.txt"),
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
gList = c("CD3D",
          "CD4","LTB","IL7R",  ##CD4 
          "CD8A","GZMB","PRF1", #CD8, NK
          "CD79A",       #Plasma
          "MS4A1",       #B
          "HLA-DQA1",    #DC from Panglao & lit "HLA-DQB1"
          "FCGR3A",       #CD16 macro> mono
          "CD14",         #mono > macro literature/schelker
          "CD68",         #macro > mono 
          "C1QA",         #macro Tirosh   
          "LILRA4",       #pDC
          "CD1C","CADM1", #cDC1 vs. cDC2 both Panglao and literature
          "S100A12",      #neutro Danaher&LM
          "CPA3",         #mast
          "MKI67")        #proliferation,Tgd    
# Get (arithmetic) average of markers: non-log scale
tmp_averages <- AverageExpression(target_seu,
        assays="RNA",
        features=gList)$RNA
tmp_averages <- t(tmp_averages)
tmp_zscores <- scale(tmp_averages)
mks_clust_avg <- cbind(rownames(tmp_zscores),tmp_zscores)

## Step 3.3: export output from enrichment and marker expression
Ftab_out <- cbind(mks_clust_avg,Ftab_candidate)
fout <- c("data/GGO_cd45pos_pred.txt")
write.table(Ftab_out, fout, sep="\t",
          quote = F, row.names = F)

saveRDS(cd45pos,file=c("GGO_cd45pos_seu.rds"))

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

## Heatmap of enrichment scores --> T cells vs. myeloid vs. B/plasma
p2 <- aheatmap(enrichedData$score_Mt,scale="col")

## Step 3.4: Manually inspect the 3.3 output and assign cell type:
## (1) if difference between 1st and 2nd highest score > 1, 
##         assign cell type associated with 1st highest score
## (2) else
##         assign based on cell type expressing canonical markers:
## Rules for assigning immune cell type based on enrichment and zscore
## NK:  T-cells or NK by enrichment & CD3D<0 & PRF1>0.5
## NKT: T-cells or NK by enrichment & CD3D>0 & PRF1>0.5
## CD4: T-cells by enrichment & sum(CD3D,LTB,IL7R) >2
## CD8: T-cells or NK by enrichment & sum(CD3D,CD8A)   
## B:  by enrichment & MS4A1>3
## plasma: by enrichment & 0<CD79A<2
## mast: by enrichment & CPA3>1
## pDC: by enrichment & LILRA4>1
## Neutrophil: either Macrophage or Monocyte or DC by enrichment, 
##             AND S100A12>1
## DC: either Macrophage or Monocyte or DC by enrichment, 
##     AND sum(CD1C,CADM1,HLA-DQA1>2)
## Mono: either Macrophage or Monocyte or DC by enrichment, 
##     AND (CD14>2)
## Macophage: either Macrophage or Monocyte or DC by enrichment, 
##     AND sum(FCGR3A,CD68,C1QA)
 
