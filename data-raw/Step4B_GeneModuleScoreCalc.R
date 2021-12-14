## ***********************************************************************
## Annotating clusters by subtypes/activities based on expression of
## mutually exclusive markers
## --> Module score of each cell = mean of z-scores of mutually exclusive 
## gene expression in the module
## --> Module score of ech cluster = mean of module scores across cells in
## the cluster
## Input:
## (1) seurat object of selected cells
## (2) Gene sets: genes x cell_type 
## Output: F4zmean of gene modules x cells
## run in 3.6.3. Seurat version 3.1.5 
## ***********************************************************************
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(NMF)
library(reshape)

source("data/source_fxns4seurat.R")

## specify input: (1) seurat object and 
## (2) genesets, including text file and the selected gene sets
## E.g. fibroblast and gene sets from Kieffer et al.
fin1 <- c("example/GGO_fibro_seu_gitTest.rds") 
fin_gset <- paste("data/","Kieffer_CAFs_agset.txt",sep="")
selGSet <-  c("detox_iCAF","IFNG_iCAF","IL_iCAF",
            "ecm_myCAF","TGFbeta_myCAF","wound_myCAF")  

## Step 1: import data
kieffer_db <- read.delim(fin_db1,sep="\t",stringsAsFactors=FALSE)
db_sel = kieffer_db %>% filter(cellType %in% selGSet)
colnames(db_sel)= c("Symbol","cellType")
tmp_subT = c("CAFS1sub")

fin1 <- c("GGO_fibro_seu.rds") 
target_seu <-readRDS(fin1)


## Step 2: gene module analysis
moduleScoreList <- GeneModuleScoreCalc(target_seu,db_sel)
cell_metaData <- moduleScoreList$cell_metadata  ##cID,histology
score_perCell <- moduleScoreList$scoreMt         ##geneSets x cells

## cluster-based module score calculation: mean of cell-based module scores
keep_gset <- which(!is.na(score_perCell[,1]))
tmp_avg <- aggregate(t(score_perCell[keep_gset,]),
        by=list(cell_metaData$cID),mean)
score_perCluster <- tmp_avg[,-1]
rownames(score_perCluster) = tmp_avg[,1] 

## Step 3: plotting and statistical comparison
## Step 3.1: Distribution of module scores in each cluster by violin plots. 
tmpdata  <- as.data.frame(t(score_perCell))    #*********[-2,])
tmpdata$cID <- cell_metaData$cID
tmpdata$hist <- cell_metaData$hist
                 
tmp4plot <- melt(tmpdata,id=c("cID","hist"))
names(tmp4plot) <- c("cID","hist","geneSet","zscore")

ggbase <- ggplot(tmp4plot, aes(x=cID, y= zscore,fill=cID))
p1 = ggbase+geom_violin()+
     geom_boxplot(width=0.1,outlier.shape=NA)+
     facet_grid(geneSet ~ .,scale="free") +   #vertical
     theme(text=element_text(size=300),
          legend.text = element_text(size = 10))+
     theme_classic()
p1

## Step 3.2: Heatmap of cluster-based module scores
p2 <- aheatmap(score_perCluster,scale='col',Rowv=NA,Colv=NA,
          border=T)
p2

## Step 3.3: statistical analysis 
## by kruskal test: module_score ~ clusterID
tmpcgs<- names(table(tmp4plot$geneSet))
pv_mod_clusters <- rep(NA,length(tmpcgs))
for (i in c(1:length(tmpcgs))){
   tmpi <- which(tmp4plot$geneSet==tmpcgs[i])
   pv_mod_clusters[i]= kruskal.test(tmp4plot$zscore[tmpi]~tmp4plot$cID[tmpi])$p.value
}
names(pv_mod_clusters) <- tmpcgs


## pair-wise comparison
selClusterPair <-  c(11,9)  #****specify pair
tmp_pair_pv <- rep(1,nrow(score_perCell))
tmpdata_sub <- tmpdata %>% filter(cID %in% selClusterPair)
id <- rownames(score_perCell)
for (i in c(1:nrow(score_perCell))){
  tmp_pair_pv[i]=kruskal.test(tmpdata_sub[,id[i]]~tmpdata_sub$cID)$p.value
}
names(tmp_pair_pv) <- rownames()





