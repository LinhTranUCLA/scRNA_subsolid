# Analysis of the single-cell transcriptome of early lung adenocarcinoma

Custom R script to analyze the 10X scRNA-seq of lung subsolid and solid nodules. The R script was built on Seurat package.

# Description
The project aims to:

1. Annotate cell clusters by immune and non-immune cell types The R script encodes a pipeline built on a two-tier approach, including (1) enrichment approach based on marker gene sets and (2) expression of canonical genes. 

2. Determine DEGs between malignant and matched normal tissue. DEG pipelines utilize single-cell MAST and pseudo-bulk edgeR algorithms. The DEG analyses incorporate patient ID as a variable in modeling.

3. Determine malignant-associated ligand-receptor interactions between cell lineages. 
