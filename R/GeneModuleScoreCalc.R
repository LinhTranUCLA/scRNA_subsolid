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

source("source_fxn4seurat.R")

## Step 0: setting up working directory
workFolder = c("scRNA_nodules/Ranalysis_Gencode34/")
setwd(workFolder)

## Get gene sets, i.e. Kieffer et al. for fibroblast
db_folder <- c("scRNA_nodules/immunedb/")
fin_db1 <- paste(db_folder,"Kieffer_CAFs_agset.txt",sep="")
kieffer_db <- read.delim(fin_db1,sep="\t",stringsAsFactors=FALSE)

## select gene sets of interest
selGSet = c("detox_iCAF","IFNG_iCAF","IL_iCAF",
            "ecm_myCAF","TGFbeta_myCAF","wound_myCAF")  
db_sel = kieffer_db %>% filter(cellType %in% selGSet)
colnames(db_sel)= c("Symbol","cellType")
tmp_subT = c("CAFS1sub")

## Step 1: load Seurat object of selected cells, e.g. AT2 cells
fin1 <- c("GGO_target_seu.rds") 
target_seu <-readRDS(fin1)

## Step 2: gene module analysis
## Step 2.1: identify shared genes between expression data and database
tmpsel <- is.element(db_sel$Symbol,rownames(target_seu))  #*****
g4hm <- db_sel[tmpsel,]

## Step 2.2: get genes exclusively expressed based on database
tmp_rep <- table(g4hm[,1])

tmp4keep <- names(tmp_rep)[tmp_rep==1]  
g4hm <- g4hm[is.element(g4hm[,1],tmp4keep),]
gList <- unique(g4hm[,1])

## Step 2.3: z-transform and gene module score calculation
## get expression
tmp2 <- t(FetchData(target_seu,vars=gList))  

## get cell information
tmp4cell=data.frame(cID = droplevels(Idents(target_seu)),
                    hist = target_seu@meta.data$histology)                 
tmp4cell$hist = factor(target_seu@meta.data$histology,
                 levels=c("Normal","GGO","Tumor"))
rownames(tmp4cell) = colnames(tmp2)

## reorder cells based on tissue type
tmpOrder <- order(tmp4cell$cID,tmp4cell$hist)
tmp4cell <- tmp4cell[tmpOrder,]
tmp2.n <- tmp2[,tmpOrder]

## Cell-based module score calculation: z-mean per pathway in each cell
npath = length(selGSet)
F4zmean <- matrix(0,npath,ncol(tmp2.n))
for (k in c(1:npath)){
    tmpg = db_sel$Symbol[which(db_sel$cellType==selGSet[k])]
    tmpi  <- which(is.element(rownames(tmp2.n),tmpg))
    F4zmean[k,]=colMeans(tmp2.n[tmpi,,drop=FALSE])
}
colnames(F4zmean)=colnames(tmp2.n)
rownames(F4zmean)=selGSet

## cluster-based module score calculation: mean of cell-based module scores
keep_gset <- which(!is.na(F4zmean[,1]))
tmp_avg <- aggregate(t(F4zmean[keep_gset,]),
        by=list(tmp4cell$cID),mean)
score_perCluster <- tmp_avg[,-1]
rownames(score_perCluster) = tmp_avg[,1] 

## Step 3: plotting and statistical comparison
## Step 3.1: Distribution of module scores in each cluster by violin plots. 
tmpdata  <- as.data.frame(t(F4zmean))    #*********[-2,])
tmpdata$cID <- tmp4cell$cID
tmpdata$hist <- tmp4cell$hist
                 
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
selClusterPair <-  c(11,8)  #****specify pair
tmp_pair_pv <- rep(1,nrow(F4zmean))
tmpdata_sub <- tmpdata %>% filter(cID %in% selClusterPair)
id <- rownames(F4zmean)
for (i in c(1:nrow(F4zmean))){
  tmp_pair_pv[i]=kruskal.test(tmpdata_sub[,id[i]]~tmpdata_sub$cID)$p.value
}
names(tmp_pair_pv) <- rownames(F4zmean)





