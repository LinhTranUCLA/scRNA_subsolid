## *********************************************************************
## Import data (filtered_feature_bc_matrix) from cellranger
## Input is tab delimited text file with four columns: 
## 1st column lists path to "filtered_feature_bc_matrix",
## 2nd column lists sample IDs, 
## 3rd lists abbreviation of sampleIDs (i.e. included in cell barcodes), 
## 4th column is batch information
## Example: "fullPath/Case1_Normal/outs/filtered_feature_bc_matrix/" "Case1_Normal" "C1N" "batch1"
## Note: each input file per batch
## *********************************************************************
importCellRangerCntMts <- function(
   CellRangeInputDirs){
 
   inputList <- as.matrix(read.delim(CellRangeInputDirs,sep="\t",
      header=F,stringsAsFactors=FALSE))
    
   allPaths = inputList[,1]
   sampID_vec = inputList[,2]
   sampID_short = inputList[,3]
   batchInfo    = inputList[,4]

   ## import 10X data and create seurat objects
   nsamples = length(allPaths)
   tmp1 <- Read10X(data.dir = allPaths[1])
   tmp1_seu <- CreateSeuratObject(counts = tmp1,
                     project = sampID_vec[1])

   allDataRaw = list()
   for (sample.Idx in c(2:nsamples)){
      tmp2 <- Read10X(data.dir = allPaths[sample.Idx])
      tmp2_seu <- CreateSeuratObject(counts = tmp2,
                     project = sampID_vec[sample.Idx])
      allDataRaw[[(sample.Idx-1)]] <- tmp2_seu  
   }

   ## merge samples from same batch together
   batch_seu <- merge(x=tmp1_seu,y=unlist(allDataRaw),
              add.cell.ids=sampID_short,
              project=batchInfo[1])
   return(batch_seu)
}

## *********************************************************************
## Perform enrichment analysis to determine if cluster markers
## is enriched by cell type markers
## Input: (1) target_seu: seurat object with cluster id, 
##        (2) mks_tab: positive cluster markers 
##             (i.e. output of Seurat FindAllMarkers), 
##        (3) db_sel: cell markers in data frame of two variables: 
##              "Symbol" and "cellType"
## Output: Matrix of enrichment scores: clusterID x cell_types
## *******************************************************************

enrichScoreCalc_mt <- function(target_seu,mks_tab,db_sel){

   isDetected <- is.element(db_sel$Symbol,rownames(target_seu))
   db_sel <- db_sel[isDetected,]
   clusID_num = table(mks_tab$cluster)
   cellTp_num = table(db_sel$cellType)
   ngtotal = nrow(target_seu)
   Ftab <- matrix(NA,length(clusID_num),length(cellTp_num))
   Fcom <- NULL
   for (i in c(1:nrow(Ftab))){
      tmpm_c <- mks_tab$gene[mks_tab$cluster==names(clusID_num)[i]]
      for (j in c(1:ncol(Ftab))){
         tmpm_t <- db_sel$Symbol[db_sel$cellType==names(cellTp_num)[j]]
         Ftab[i,j] = length(intersect(tmpm_c,tmpm_t))*ngtotal/
                 clusID_num[i]/cellTp_num[j]  
         tmpo <- intersect(tmpm_c,tmpm_t)
         Fcom = rbind(Fcom,cbind(tmpo,
                 rep(names(cellTp_num)[j],length(tmpo)),
                 rep(names(clusID_num)[i],length(tmpo))))       
       }
    }

    rownames(Ftab) =names(clusID_num)
    colnames(Ftab) =names(cellTp_num)
    colnames(Fcom)= c("symbol","nM_id","clusterid")
    list(score_Mt=Ftab,overlapMks=Fcom)
}

## ***************************************************************
## Get sum of raw reads of all cells per sample as bulk expression
## Input: expression is genes x cells, and 
## information of sample cells are associated with
## ***************************************************************
getColGroupSum <- function(x,groupid){
   ngr = length(groupid)
   output <- matrix(0,nrow(x),ngr)
   for (i in c(1:ngr)){
      tmpsel = which(colnames(x)==groupid[i])
      output[,i]=rowSums(x[,tmpsel])
   }
   colnames(output) = groupid
   rownames(output) = rownames(x)
   return(output)
}

## ********************************************************************
## Calculate (genometric) average (log) expression per each cell type in 
## each sample
## Output is a list whose element is a matrix of gene expression x samples
## for each cell type
## *********************************************************************
getExp_CellType_PerSample <- function(target_seu,selCellTypes_info){

   ## get cell type of interest
   list_cellType <- selCellTypes_info[!is.na(selCellTypes_info[,2]),2]
   list_cellType <- unique(list_cellType)
   
   ## initiate list of expr avg per cell type.
   ncellTypes <- length(list_cellType)
   Exp_avg <- vector("list",ncellTypes)
   names(Exp_avg) <- list_cellType

   # Create average matrix per sample per each cell type
   for (k in c(1:ncellTypes)){
      tmpidx <- which(selCellTypes_info[,2]==list_cellType[k])
      tmpselClus <- selCellTypes_info[tmpidx,1]
      subData_seu <- subset(target_seu,idents=tmpselClus)
      tmp_data <- as.matrix(GetAssayData(subData_seu, 
             assay= "RNA",slot = "data"))
      tmp_avg <- aggregate(t(tmp_data),
             by=list(subData_seu@meta.data$orig.ident),mean)
      tmp_avg2 <- tmp_avg[,-1]
      rownames(tmp_avg2) = tmp_avg[,1]
      Exp_avg[[k]] <- t(tmp_avg2)
   }
   return(Exp_avg)
}

## ******************************************************************
## Generate color vector for clusters
## ******************************************************************
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

## ************************************************************************
## Calculate Gene Module Score per cell
## Gene sets will be filtered to obtain the mutual exclusive markers
## before calulating scores.
## --> Module score of each cell = mean of z-scores of mutually exclusive 
## gene expression in the module
## --> Module score of ech cluster = mean of module scores across cells in
## the cluster
## Input files
## (1) seurat object of selected cells with "histology" in meta-data
## (2) data frame of (a)Symmbol and (b) cellType
## Output: list of (1)meta-data of cells and (2) F4zmean of gene modules x cells
## *************************************************************************
GeneModuleScoreCalc <- function(target_seu,db_sel){
   
   selGSet <- names(table(db_sel$cellType))

   ## Step 1: identify shared genes between expression data and database
   tmpsel <- is.element(db_sel$Symbol,rownames(target_seu))  #*****
   g4hm <- db_sel[tmpsel,]

   ## Step 2: get genes exclusively expressed based on database
   tmp_rep <- table(g4hm[,1])
   tmp4keep <- names(tmp_rep)[tmp_rep==1]  
   g4hm <- g4hm[is.element(g4hm[,1],tmp4keep),]
   gList <- unique(g4hm[,1])

   ## Step 3: z-transform and gene module score calculation
   ## get expression
   tmp2 <- t(FetchData(target_seu,vars=gList))  

   ## get cell information
   tmp4cell=data.frame(cID = droplevels(Idents(target_seu)),
                    hist = target_seu@meta.data$histology)                 
   ##tmp4cell$hist = factor(target_seu@meta.data$histology,
   ##              levels=c("Normal","GGO","Tumor"))
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
   
   ModuleScoresList <- list(cell_metadata=tmp4cell,scoreMt=F4zmean)
   return(ModuleScoresList)
}
   
