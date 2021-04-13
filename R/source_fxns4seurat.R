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



