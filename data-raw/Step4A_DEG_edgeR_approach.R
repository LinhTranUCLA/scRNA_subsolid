## *************************************************************************
## DEG between normal and tumor --> using pseudo-bulk approach with edgeR
## model: express ~ patientID + histology
## Input is seurat object with information of sampleID (orig.ident) should 
## have format as PatientID_TissueType. "Normal" is used to indicate the 
## control group
## run in 3.6.3. Seurat version 3.1.5 
## *************************************************************************
library(Seurat)
library(dplyr)
library(edgeR)
library(scater)

source("data/source_fxns4seurat.R")

## Step 1: load Seurat object of the selected cells, e.g. AT2 cells
fin1 <- c("GGO_AT2_seu.rds") 
target_seu <-readRDS(fin1)

## Step 2: create pseudo bulk expression
## Get raw count
tmp_cnt <- as.matrix(GetAssayData(target_seu, 
             assay= "RNA",slot = "counts"))
colnames(tmp_cnt) <- cell_clin$orig.ident   #****change to sampleID

## remove undetected genes: detected in <1% of cells
tmp_detect <- rowSums(sign(tmp_cnt))
#summary(tmp_detect)
keptR <- which(tmp_detect>round(0.01*ncol(tmp_cnt)))  #***detected in 1% cells
tmp_cnt = tmp_cnt[keptR,]

# create pseudo-bulk data
tmp_bulk <- getColGroupSum(x=tmp_cnt,groupid=unique(colnames(tmp_cnt)))

# extract sample information based on sample ID
bulk_group = t(sapply(colnames(tmp_bulk),function(x) unlist(strsplit(x,"_")),
          simplify=T))
colnames(bulk_group) = c("Case","Histo")
isMalig <- ifelse(bulk_group[,2]=="Normal",0,1)
bulk_group <- cbind(bulk_group,isMalig)

## Step 3: DEG by edgeR
## Step 3.1: create object
d1 = DGEList(genes=rownames(tmp_bulk),counts=tmp_bulk,samples=bulk_group)

## Step 3.2: preprocessing
## Step 3.2.1: filter "samples" with low library 
## (i.e. samples with small number of cells)
bad_samples <- isOutlier(d1$samples$lib.size,log=T,type="lower")
bad_samples

## Step 3.2.1: remove genes with (1) low rc and (2) only expressed in single sample
## by edgeR "filterByExpr". 
thres_cnt <- ceiling(0.5*min(table(colnames(tmp_cnt))))
tmpTest = filterByExpr(d1, group = d1$samples$isMalig,min.count=thres_cnt)
keepG = tmpTest  #****

## Step 3.2.1: Normalization by M-value method
d1 <- d1[keepG,]
d1 <- calcNormFactors(d1)
d1$samples

## Step 3.3: statistical model by edgeR
## Step 3.3.1: create model for regression
design <- model.matrix(~factor(Case) + factor(isMalig), d1$samples)
design

## Step 3.3.2: estimate the negative binomial (NB) dispersions
d1 <- estimateDisp(d1, design)  
summary(d1$trended.dispersion)
plotBCV(d1)

## Step 3.3.3: fit a GLM with Quasi-likelihood Tests
fit <- glmQLFit(d1, design, robust=TRUE)
summary(fit$var.prior)

## Step 3.4: statitical test and export 
res <- glmQLFTest(fit, coef=ncol(design))  #the last factor in design = isMalig
summary(decideTests(res))

# extract all DEGs
p_thres=0.05
lfc_thres=log2(1.5)

tmp_degs = topTags(res,n=nrow(res),adjust.method="BH")

fout <- paste("DEG_byEdgeR_p0.05_fc1.5.txt")
write.table(tmp_degs,fout,sep="\t",quote=F,row.names=F)
