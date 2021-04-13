## *************************************************************************
## Example of DEG between normal and tumor by using MAST approach
## express ~ patientID + histology
## Input:
## (1) seurat object:
##    @meta.data$orig.iden: sampleID as: PatID_TissueLesion
## 	@meta.data$binaryHis: tissue types, reference/control = "Normal"
##	@meta.data$PatID: patient ID (not sample ID) due to 
##		multiple samples/patient 
## (2) txt file of two columns: 
## 	col 1: cluster ID 
##	col 2: associated cell types. Use "NA" to indicate sample specific
##          clusters --> not included in the analysis 
## Output: text files list DEGs genes and related statistical values
## File name = paste("KvN_select",cellType,"_mast.txt",se[="")
## run in 3.6.3. Seurat version 3.1.5 
## *************************************************************************
library(MAST)
library(Seurat)
library(dplyr)    
library(limma)
library(reshape2)
library(data.table)
library(knitr)

## Step 0: setting up working directory
workFolder = c("scRNA_nodules/Ranalysis_Gencode34/")
setwd(workFolder)

source("source_fxn4seurat.R")

## Step 1: load Seurat object of selected cells, e.g. cd45neg cells
fin1 <- c("GGO_cd45neg_seu.rds") 
cd45neg <-readRDS(fin1)  

fin2 <- c("GGO_cd45neg_selClusters.txt")  #***w/disease status
cellAnnotInfo_cd45neg <- read.delim(fin2,sep="\t",header=T,
        stringsAsFactors=FALSE)

## Step 2: DEGs by MAST - loop
## parameter for MAST
bin_no <- 20 #*****for thresholdSCRNACountMatrix
bin_mem <- 30 #****for thresholdSCRNACountMatrix
thrs_freq1 <- 0.20 #*****filter out poorly expressed genes 

## create inputs
target_seu <- cd45neg
DefaultAssay(target_seu) <- 'RNA'
tmp_cellAnnotInfo <- cellAnnotInfo_cd45neg

# loop calculate Cancer vs N DEGs for each cell type
tmp <- tmp_cellAnnotInfo[,2]  #list celltype of selected clusters 
cellTypeID <- tmp[!duplicated(tmp) & !is.na(tmp)]
ncellTypes <- length(cellTypeID)
for (k in c(1:ncellTypes)){
   tmpselType = cellTypeID[k]
   tmpi <- which(tmp_cellAnnotInfo[,2]==tmpselType)
   tmpselClus <- tmp_cellAnnotInfo[tmpi,1]
   tmpselCell <- which(is.element(Idents(target_seu),tmpselClus))
   subData_test  <- subset(target_seu,cells=tmpselCell)
   ##subData_test[["binaryHis"]] <- relevel(subData_test@meta.data$binaryHis,
   ##                            ref=c("Normal"))   #******

   tmp_clin <- FetchData(subData_test,vars=c("PatID","binaryHis"))  #***
   tmp_clin$PatID <- factor(tmp_clin$PatID)
   tmp_clin$histology <- relevel(tmp_clin$binaryHis,
             ref=c("Normal"))
   tmp_cnt <- as.matrix(GetAssayData(subData_test,
             assay= "RNA",slot = "counts"))

   ## remove poorly expressed genes: < thrs_freq1
   gene_sparsity <- rowSums(tmp_cnt!=0)
   filter_by_sparse_genes <- gene_sparsity >=thrs_freq1*ncol(tmp_cnt)  #*****

   tmp_cnt <- tmp_cnt[filter_by_sparse_genes,]
   tmp_fea <- data.frame(GSymbol=rownames(tmp_cnt))

   tmp_norm <- as.matrix(NormalizeData(tmp_cnt,
              normalization.method = "LogNormalize", scale.factor = 1e6))
   sca_all <- FromMatrix(tmp_norm,
         cData=as.data.frame(tmp_clin),
         fData=tmp_fea, class = "SingleCellAssay")

   thres <- thresholdSCRNACountMatrix(assay(sca_all), 
         nbins = bin_no , min_per_bin = bin_mem)
   assays(sca_all) <- list(thresh=thres$counts_threshold, 
         tpm=assay(sca_all))  #store bin ($thres) which matrix element belong

   zlmCond <- zlm(~ binaryHis + PatID, sca_all)  #*******
   selcontrast <- colnames(zlmCond@coefC)[2]
   summaryCond <- summary(zlmCond, doLRT=selcontrast)
   summaryDt_cond <- summaryCond$datatable
   fcHurdle_cond <- merge(summaryDt_cond[contrast==selcontrast & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
      summaryDt_cond[contrast==selcontrast & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], 
      by='primerid') #logFC coefficients
   fcHurdle_cond[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] ## add fdr column 
   tmpo <- order(-abs(fcHurdle_cond$coef))

   # store the DEG results in "SubTypes" folder
   tmpfout = paste("KvN_selected",tmpselType,"_mast.txt",sep="")
   write.table(fcHurdle_cond[tmpo,],tmpfout,sep="\t",row.names=F,quote=F)
}



