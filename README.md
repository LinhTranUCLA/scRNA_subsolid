# Analysis of the single-cell transcriptome of early lung adenocarcinoma

The repository included a custom R script to analyze the 10X scRNA-seq of lung subsolid and solid nodules. The single-cell objects were built on the Seurat package.

# Description
The project aims to:

1. Annotate cell clusters for immune and non-immune cell types. The R script encodes a pipeline built on a two-tier approach, including (1) enrichment approach based on marker gene sets and (2) expression of canonical genes. 

2. Determine DEGs between malignant and matched normal tissue. DEG pipelines utilize single-cell MAST and pseudo-bulk edgeR algorithms. The DEG analyses incorporate patient ID as a variable in modeling.

3. Determine malignant-associated ligand-receptor interactions between cell lineages. 


# Data
All cellranger count matrices and clinical data are available at rds files are available at HTAN, HUMAN TUMOR ATLAS NETWORK,
https://humantumoratlas.org/HTA3/

The rds files for immune and non-immune cells are available at:
https://drive.google.com/drive/folders/1OYIltk5TeoNTKEfGvkpNg_1Z_-XeterL?usp=sharing

# Note
The step-by-step pipeline was illustrated by scripts located in the data-raw directory. The inputs of the examples will be found in the example directory.
