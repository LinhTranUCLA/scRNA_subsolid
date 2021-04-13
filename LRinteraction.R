## ***********************************************************************
## Identify the interactions among cell types
## Ligand-receptor information is obtained through iTALK package
## Input:
## (1) seurat object:
## 	@meta.data$binaryHis: tissue types/groups, reference = "Normal"
##	@meta.data$PatID: patient ID (not sample ID) due to 
##		multiple samples/patient 
## (2) txt file of two columns: 
## 	col 1: cluster ID 
##	col 2: associated cell types. Use "NA" to indicate sample specific
##          clusters --> not included in the analysis 
## (3) txt files of DEGs in each cell type: output from 
##    "example_DEG_MAST_multipleCellType.R"
## run in 3.6.3. Seurat version 3.1.5 
## ************************************************************************  
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(reshape2)  
library(NMF)
library(iTALK)     

# create index of lower triangle (include diagonal) of n x n matrix
lowerTriInd_wDiag <- function(n) {
  z <- sequence(n)
  cbind(
    row = unlist(lapply(1:n, function(x) x:n), use.names = FALSE),
    col = c(rep(z[-length(z)], times = rev(tail(z,-1))),n))
}

# create index of lower triangle (no diagonal) of n x n matrix
lowerTriInd_noDiag <- function(n) {
  z <- sequence(n)
  cbind(
    row = unlist(lapply(2:n, function(x) x:n), use.names = FALSE),
    col = rep(z[-length(z)], times = rev(tail(z,-1))-1))
}


# get statistical infor from kruskal test
getKruskalResults <- function(x,pheno){
   x_mean <- aggregate(x,by=list(pheno),mean)[,-1]
   x_pv <- kruskal.test(x~pheno)$p.value
   c(x_mean,x_pv)
}

getZscore_mean <- function(x,pheno){
   mean_x = mean(x[pheno=="Normal"],na.rm=T)
   sd_x   = sd(x[pheno=="Normal"],na.rm=T)
   if (sd_x!=0){
      zscore = (x-mean_x)/sd_x
      zscore_mean = mean(zscore[pheno!="Normal"])
   }else{
      zscore_mean = 0
   }
   return(zscore_mean)
}

## Step 0: setting up working directory
workFolder = c("scRNA_nodules/Ranalysis_Gencode34/")
setwd(workFolder)

## Step 1: load Seurat object of selected cells, e.g. cd45neg cells
fin1 <- c("GGO_cd45neg_seu.rds") 
cd45neg <-readRDS(fin1)  

fin2 <- c("GGO_cd45neg_selClusters.txt")  #***w/disease status
cellAnnotInfo_cd45neg <- read.delim(fin2,sep="\t",header=T,
        stringsAsFactors=FALSE)

fin3 <- c("GGO_cd45pos_seu.rds") 
cd45neg <-readRDS(fin3)  

fin4 <- c("GGO_cd45pos_selClusters.txt")  #***w/disease status
cellAnnotInfo_cd45pos <- read.delim(fin4,sep="\t",header=T,
        stringsAsFactors=FALSE)

## prefix and suffix of DEG files: prefix_cellType_suffix, 
## prefix includes full path
degFile_prefix <- c("KvN_select")
degFile_suffix <- c("_mast.txt")

## setup input
DefaultAssay(cd45neg) <- 'RNA'
DefaultAssay(cd45pos) <- 'RNA'

tmpClust = Idents(cd45neg)
cd45neg <- AddMetaData(cd45neg,
           metadata=tmpClust,col.name='cd45neg_clustID')

tmpClust = Idents(cd45pos)
cd45pos <- AddMetaData(cd45pos,
           metadata=tmpClust,col.name='cd45pos_clustID')

## Step 2: average (log) expression per each cell type per sample
Exp_cd45neg <- getExp_CellType_PerSample(cd45neg,cellAnnotInfo_cd45neg)
Exp_cd45pos <- getExp_CellType_PerSample(cd45pos,cellAnnotInfo_cd45pos)
Exp_all <- c(Exp_cd45neg,Exp_cd45pos)

## check cells detected in all samples
nsamples <- length(table(cd45neg@meta.data$orig.iden))
mm <- sapply(Exp_all,function(x) ncol(x))
kept_cellType <- mm==nsamples

Exp_all <- Exp_all[kept_cellType]

## Step 3: calculate interaction scores of between pair of cell type
## per each sample
## Step 3.1: extract sample information from sampleID
tmpID <- colnames(Exp_all[[1]])
tmp_his = t(sapply(tmpID,function(x) unlist(strsplit(x,"_")),
          simplify=T))

tmpidx = grep("GGO",tmp_his[,2])
tmp_his[tmpidx,2]="GGO"
binHis = ifelse(tmp_his[,2]!="Normal","GGO_Tu",tmp_his[,2])

SampleInfo <- data.frame(SampleID = tmpID, PatID=tmp_his[,1],
               histology=tmp_his[,2],
               binaryHis=binHis)
SampleInfo$histology <- factor(SampleInfo$histology,
              levels=c("Normal","GGO","Tumor"))
SampleInfo$binaryHis <- factor(SampleInfo$binaryHis,
              levels=c("Normal","GGO_Tu"))
 
## Step 3.2: load LR pair r=from iTALK package
database <- iTALK:::database

## Step 3.3: LR pairs with available expression 
allGenes = data.frame(symbol=rownames(Exp_all[[1]]),
                      index=c(1:nrow(Exp_all[[1]])))

## select type of LR interaction: "checkpoint","cytokine,"growthfactor","other"
list_LRtype <- c("checkpoint","cytokine","other")

## determining possible cell type pairs
list_cellType <- names(Exp_all)
ncellTypes <- length(list_cellType)

index_inter <- lowerTriInd_wDiag(ncellTypes)   #2cols matrix
cID_inter <- cbind(list_cellType[index_inter[,1]],
                   list_cellType[index_inter[,2]])
colnames(cID_inter) <- c("fromCell","toCell")
name_inter <- apply(cID_inter,1,
      function(x) paste0(x,collapse="_"))
no_inter <- ncellTypes*(ncellTypes+1)/2

## Step 3.3: calculate LR interaction scores for pairs of cell types
## in both directions: A-to-B and B-to-A
## including (1) performng statistical tests comparing scores between disease vs. Normal,
## (2) log(FC) beween disease vs. Normal of ligands and receptors in related cell types
## (3) calculating correlation coefficient between ligand in A and receptor in B, and vice versa
Ftab_sum <- NULL  #data frame from stacked data
nsamples <- ncol(Exp_all[[1]])

for (pp in c(1:length(list_LRtype))){
   type_LR = list_LRtype[pp]

   #mapping database to expression
   ligand_ind  <-which(database$Ligand.ApprovedSymbol %in% allGenes$symbol &
                    database$Classification %in% type_LR)
   receptor_ind<-which(database$Receptor.ApprovedSymbol %in% allGenes$symbol &
                   database$Classification %in% type_LR)
   ind<-intersect(ligand_ind,receptor_ind)  #both lignad and receptor covered

   map_table <- database[ind,c('Pair.Name','Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
      left_join(allGenes,by=c('Ligand.ApprovedSymbol'='symbol')) %>%
      dplyr::rename(ligand_index=index) %>%
      left_join(allGenes,by=c('Receptor.ApprovedSymbol'='symbol')) %>%
      dplyr::rename(receptor_index=index)

   # calculate interaction of all cell type pairs
   for (k in c(1:no_inter)){
       indA <- index_inter[k,1]
       indB <- index_inter[k,2]
       dataA <- Exp_all[[indA]]
       dataB <- Exp_all[[indB]]

       # load DEGs between GT vs. N of the selected cell types
       fin_A <- paste(degFile_prefix,cID_inter[k,1],degFile_suffix,sep="")
       deg_A <- read.delim(fin_A,sep="\t",header=T,stringsAsFactors=FALSE)

       fin_B <- paste(degFile_prefix,cID_inter[k,2],degFile_suffix,sep="")
       deg_B <- read.delim(fin_B,sep="\t",header=T,stringsAsFactors=FALSE)
       
       # score matrices in which each row is interaction between a pair LR per pair A-B cells       
       tmp_A2B <- dataA[map_table$ligand_index,,drop=FALSE]*
                  dataB[map_table$receptor_index,,drop=FALSE]
       rownames(tmp_A2B)=map_table$Pair.Name 

       # statistical tests (kruskal test) comparing scores between disease and Normal       
       stat_A2B <- t(apply(tmp_A2B,1,getKruskalResults,
              pheno=SampleInfo$binaryHis))
       colnames(stat_A2B) = c("avg_nor","avg_GT","pv")
       stat_A2B <- as.data.frame(stat_A2B)
       stat_A2B$cell_from = as.character(rep(cID_inter[k,1],nrow(tmp_A2B)))  #****
       stat_A2B$cell_to = as.character(rep(cID_inter[k,2],nrow(tmp_A2B)))    #****
       stat_A2B$Classification = as.character(type_LR,nrow(tmp_A2B))
       stat_A2B$Pair.Name = map_table$Pair.Name
       stat_A2B$Ligand.ApprovedSymbol=map_table$Ligand.ApprovedSymbol
       stat_A2B$Receptor.ApprovedSymbol=map_table$Receptor.ApprovedSymbol

       lig_rat <-rowMeans(dataA[map_table$ligand_index,SampleInfo$binaryHis=="GGO_Tu",drop=FALSE],na.rm=T)-
              rowMeans(dataA[map_table$ligand_index,SampleInfo$binaryHis=="Normal",drop=FALSE],na.rm=T)
       rep_rat <-rowMeans(dataB[map_table$receptor_index,SampleInfo$binaryHis=="GGO_Tu",drop=FALSE],na.rm=T)-
              rowMeans(dataB[map_table$receptor_index,SampleInfo$binaryHis=="Normal",drop=FALSE],na.rm=T)
       stat_A2B$cell_from_mean_exprs <- lig_rat  ##lig_zscore ## lig_rat   
       stat_A2B$cell_to_mean_exprs <- rep_rat   ##rep_zscore   ##rep_rat
       stat_A2B$cor <- apply(cbind(dataA[map_table$ligand_index,,drop=FALSE],
               dataB[map_table$receptor_index,,drop=FALSE]),1,function(x)
               cor(x[1:nsamples],x[(nsamples+1):(2*nsamples)],use = "pairwise.complete.obs",method = c("spearman")))
       tmp_degAB <- stat_A2B[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
            left_join(deg_A[,c('primerid','fdr')],by=c('Ligand.ApprovedSymbol'='primerid')) %>%
            dplyr::rename(ligand_fdr=fdr) %>%
            left_join(deg_B[,c('primerid','fdr')],by=c('Receptor.ApprovedSymbol'='primerid')) %>%
            dplyr::rename(receptor_fdr=fdr)
       stat_A2B$ligand_fdr <- tmp_degAB$ligand_fdr
       stat_A2B$receptor_fdr <- tmp_degAB$receptor_fdr
            
       tmp_B2A <- dataB[map_table$ligand_index,,drop=FALSE]*
                  dataA[map_table$receptor_index,,drop=FALSE]
       rownames(tmp_B2A)=map_table$Pair.Name

       stat_B2A <- t(apply(tmp_B2A,1,getKruskalResults,
            pheno=SampleInfo$binaryHis))
       colnames(stat_B2A) = c("avg_nor","avg_GT","pv")
       stat_B2A <- as.data.frame(stat_B2A)
       stat_B2A$cell_from = as.character(rep(cID_inter[k,2],nrow(tmp_B2A)))  #*****
       stat_B2A$cell_to = as.character(rep(cID_inter[k,1],nrow(tmp_B2A)))   #****
       stat_B2A$Classification = as.character(type_LR,nrow(tmp_B2A))
       stat_B2A$Pair.Name = map_table$Pair.Name
       stat_B2A$Ligand.ApprovedSymbol=map_table$Ligand.ApprovedSymbol
       stat_B2A$Receptor.ApprovedSymbol=map_table$Receptor.ApprovedSymbol

       lig_rat <-rowMeans(dataB[map_table$ligand_index,SampleInfo$binaryHis=="GGO_Tu",drop=FALSE])-
             rowMeans(dataB[map_table$ligand_index,SampleInfo$binaryHis=="Normal",drop=FALSE])
       rep_rat <-rowMeans(dataA[map_table$receptor_index,SampleInfo$binaryHis=="GGO_Tu",drop=FALSE])-
              rowMeans(dataA[map_table$receptor_index,SampleInfo$binaryHis=="Normal",drop=FALSE])
       stat_B2A$cell_from_mean_exprs <- lig_rat   ##lig_rat   
       stat_B2A$cell_to_mean_exprs <- rep_rat      ###rep_rat
       stat_B2A$cor <- apply(cbind(dataB[map_table$ligand_index,,drop=FALSE],
               dataA[map_table$receptor_index,,drop=FALSE]),1,function(x)
               cor(x[1:nsamples],x[(nsamples+1):(2*nsamples)],use = "pairwise.complete.obs",method = c("spearman")))
       tmp_degAB <- stat_B2A[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%   #fixed bug on 10/2020
            left_join(deg_B[,c('primerid','fdr')],by=c('Ligand.ApprovedSymbol'='primerid')) %>%
            dplyr::rename(ligand_fdr=fdr) %>%
            left_join(deg_A[,c('primerid','fdr')],by=c('Receptor.ApprovedSymbol'='primerid')) %>%
            dplyr::rename(receptor_fdr=fdr)
       stat_B2A$ligand_fdr <- tmp_degAB$ligand_fdr
       stat_B2A$receptor_fdr <- tmp_degAB$receptor_fdr
      
       Ftab_sum <- rbind(Ftab_sum,stat_A2B,stat_B2A)
   }   
}

saveRDS(Ftab_sum,"LRInteract_allEdges.rds")

## Step 4: filter out insignificant interactions: 
## e.g. (1) LR_score-based pv < 0.05, (2) L and R are not DEGs in releated cell types
Ftab_sum3 <- Ftab_sum %>% filter(pv<0.05 & 
     !(is.na(ligand_fdr) & is.na(receptor_fdr)) &
      (avg_nor!=0 & avg_GT!=0)) # & cell_to != cell_from)
Ftab_sum3$ligand_fdr =ifelse(is.na(Ftab_sum3$ligand_fdr),1,Ftab_sum3$ligand_fdr)
Ftab_sum3$receptor_fdr =ifelse(is.na(Ftab_sum3$receptor_fdr),1,Ftab_sum3$receptor_fdr)
Ftab_sum3$tmpdeg <- (Ftab_sum3$avg_GT-Ftab_sum3$avg_nor)/Ftab_sum3$avg_nor   #relative change in LT scores between groups 



## Step 5: Visualize networks: using iTALK package
## Step 5.1: setup color
colCode_all = structure(c("#A52A2A","#FFA07A","#4169E1","#8A2BE2","#40E0D0",
             "#FFA500","#F4A460","#D2691E","#D2B48C","#008000", 
             "#32CD32","#EE82EE","#DB7093","#00BFFF","#9ACD32","#FF6347"),
           names=names(table(Ftab_sum3$cell_from)))

colnames(Ftab_sum3)
#selVars <-c(8,9,10,11,4,5,6)  #****
selVars <-c(8,9,15,4,5,6,7)  #****

colCode_all2 <- colCode_all
names(colCode_all2)[9]="Neu."
names(colCode_all2)[12]="Plas."
names(colCode_all2)[2]="B"
names(colCode_all2)[5]="EC"
names(colCode_all2)[6]="CAF"
names(colCode_all2)[15]="gdT"

## Step 5.2: significant interactions among non-immnune components, i.e. AT2,fibro and endo
tmpi <- Ftab_sum3$pv < 0.05 &
        Ftab_sum3$cor > 0.3 &
        Ftab_sum3$tmpdeg > 1 &   ###2-folde change of score
        Ftab_sum3$cell_from != Ftab_sum3$cell_to &
        (Ftab_sum3$cell_from %in% c("AT2","endo","fibro")  &  
         Ftab_sum3$cell_to %in% c("AT2","endo","fibro")) 
res_cat <- Ftab_sum3[which(tmpi==TRUE),selVars]
tmpid_combo <- apply(res_cat[,c(4,5,7)],1,function(x) paste0(x,collapse="_"))
res_cat <- res_cat[!duplicated(tmpid_combo),]  ##remove duplicated interactions in databases

## interaction networks
NetView(res_cat,col=colCode_all2,
       vertex.label.cex=1,arrow.width=1,edge.max.width=5,
       label=TRUE,vertex.size=50)

# count interactions and heatmap: ligands (rows) to receptors (cols)
tmp_sum <- table(res_cat$cell_from,res_cat$cell_to) 
aheatmap(tmp_sum,scale='none')
