## *************************************************************
## DEG between normal and tumor --> using MAST approach
## express ~ patientID + tissue type
## Input is seurat object with information of :
## (1) patientID: stored in xx@meta.data$PatID. May have >2 samples/patient
## (2) tissue type: stored in xx@meta.data$binaryHis with reference
##         level indicated as "Normal"
## run in 3.6.3. Seurat version 3.1.5 
## *************************************************************
suppressPackageStartupMessages({
    library(MAST)
    library(Seurat)
    library(dplyr)    
    library(ggplot2)
    library(limma)
    library(reshape2)
    library(data.table)
    library(knitr)
    library(stringr)
})

## Step 0: setting up working directory
workFolder = c("scRNA_nodules/Ranalysis_Gencode34/")
setwd(workFolder)

## Step 1: load Seurat object of selected cells, e.g. AT2 cells
fin1 <- c("GGO_target_seu.rds") 
target_seu <-readRDS(fin1)

## Step 2: create sce object for MAST
## Step 2.1: filter out poor genes and normalize data 
tmp_cnt <- as.matrix(GetAssayData(target_seu, 
             assay= "RNA",slot = "counts"))

## remove poorly expressed genes: < 2% of cells
gene_sparsity <- rowSums(tmp_cnt!=0)
summary(gene_sparsity)
filter_by_sparse_genes <- gene_sparsity >=0.02*ncol(tmp_cnt)  #*****
tmp_cnt <- tmp_cnt[filter_by_sparse_genes,]

## normalize count (seurat)
tmp_norm <- as.matrix(NormalizeData(tmp_cnt,
              normalization.method = "LogNormalize", scale.factor = 1e6))

## create SingleCellAssay (not SCE) whose wellKey in colData is cell barcode
sca_all <- FromMatrix(tmp_norm,
    cData=as.data.frame(cell_clin,stringsAsFactors = FALSE),
    fData=tmp_fea, class = "SingleCellAssay")

colnames(rowData(sca_all))
colnames(colData(sca_all))

# Step 2.2 : calculate library size of each cell
cdr2 <- colSums(assay(sca_all)>0)
colData(sca_all)$cngeneson <- scale(cdr2)

## Step 3: Adaptive thresholding
## setting parameters
bin_no <- 20 #*****
bin_mem <- 30 #****
freq_expressed <- 0.20 #******

thres <- thresholdSCRNACountMatrix(assay(sca_all), 
         nbins = bin_no , min_per_bin = bin_mem)
##length(thres$valleys)

## add bin data to sca object
assays(sca_all) <- list(thresh=thres$counts_threshold, 
       tpm=assay(sca_all))  #store bin ($thres) which matrix element belong
tmpfreq <- freq(sca_all)
##summary(tmpfreq)

## filter genes AGAIN 
expressed_genes <- freq(sca_all) > freq_expressed  
sca_all <- sca_all[expressed_genes,] 

## Step 4: Differential Expression using a Hurdle model
## Step 4.1: specify factor to perform regression
cond<-colData(sca_all)$binaryHis  #*******must be factor
cond<-relevel(cond,"Normal")      #*******normal is reference

table(cond)  #check refernce level

colData(sca_all)$condition<-cond
colnames(colData(sca_all))

## Step 4.2: fit the models: **** cond + patID
zlmCond <- zlm(~ condition + PatID, sca_all)

## Step 4.3: statistical evaluation by LRT test
tmp_stimulate <- names(table(cond))[2]
coef_id=paste("condition",tmp_stimulate,sep="")
summaryCond <- summary(zlmCond, doLRT=coef_id) 
summaryLmer2 <- summary(lmer.output,doLRT=coef_id)

## Step 4.4: export data
fcThres = log(1.5)

cpConds = paste("condition",levels(colData(sca_all)$condition)[2],sep="")
summaryDt_cond <- summaryCond$datatable ## extract whole table

fcHurdle_cond <- merge(summaryDt_cond[contrast==cpConds & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt_cond[contrast==cpConds & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], 
                  by='primerid') #logFC coefficients
fcHurdle_cond[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] ## add fdr column 

fout <- c("DEG_byMastLRT.txt")
write.table(fcHurdle_cond,fout,sep="\t",quote=F,row.names=F)

