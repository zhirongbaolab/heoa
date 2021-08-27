# HEOA: human Early Organogenesis Atlas

## 1. This repository contains the R code used in the following papers on human emrbyo:

A single cell transcriptome atlas of human early embryogenesis

## 2. The introduction for file structure

  1. *.html: R markdown files for generating each figure.
  2. src/: source codes for each figure (corresponding to *.html) and all other analysis (identification of cell types and identification of systemically changing genes).
  3. data/: necessary input data for scripts: for example, gene annotation, cell annotation, color scheme, contours of limb domains.
  4. list/: gene list from published data: batch effect genes, cell cycle genes, and hemoglobin genes.
  5. plot/: figures generated in the analysis that are not included in R markdown files.
  6. remap_Tan_frog: the script to map RNA-seq data of frog embryos (Tan et al. 2013) on latest frog genome.
  
## 3. An online depository for cell types and gene expression is avaiable at https://heoa.shinyapps.io/base/.

## 4. The raw and processed data can be downloaded from the NCBI Gene Expression Omnibus (GSE157329). 
