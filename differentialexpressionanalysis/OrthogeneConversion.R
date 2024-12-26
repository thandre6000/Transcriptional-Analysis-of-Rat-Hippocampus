
setwd("/Users/theoandre/Desktop/Deconvolution")    

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("orthogene")

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")

# load necessary packages

library(tidyverse)
library(xlsx)
library(orthogene)
library(Seurat)


#read in bulk DEseq data

Countdata <- readxl::read_xlsx("/Users/theoandre/Downloads/filtered_counts (1).xlsx")
colnames(Countdata)[1] <- "gene_id"
rownames(Countdata) <- Countdata$gene_id
rownames <- Countdata$gene_id
Countdata <- Countdata[,c(2:ncol(Countdata))]
rownames(Countdata) <- rownames
rownamesconverted <- orthogene::convert_orthologs(gene_df = Countdata,
                                                  gene_input = "rownames", 
                                                  gene_output = "rownames",
                                                  input_species = "rat",
                                                  output_species = "mouse",
                                                  non121_strategy = "keep_popular",
                                                 ) 







rownames(Countdata) <- rownames

#Convert count list to orthologs

colnames(Countdata)[1] <- "gene_id"
gene_df <- orthogene::convert_orthologs(gene_df = Countdata,
                                          +                                         gene_input = "gene_id", 
                                          +                                         gene_output = "rownames", 
                                          +                                         input_species = "rat",
                                          +                                         output_species = "mouse"
                                                                              ) 

