
title: "High-Res DEseq"
author: "Theo Andre"
date: "2023-07-27"
output: html_document

## DEseq cuz im losing my mind 

BiocManager::install("DESeq2")
install.packages("readxl")
library(readxl)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(biomaRt)

#determine library size
reads <- readxl::read_xlsx("/Users/tandre/Downloads/filtered_counts.csv")
sums <- as.data.frame(colSums(reads[,c(2:45)]))
View(as.data.frame(sums))

#REVERSE THESE AND CLEAN UP

sums[,2] <- rownames(sums)
sums1 <- as.data.frame(sums[rownames(sums) %in% rownames(design1) ,])
rownames(sums1) <- rownames(sums)

#check if any group has higher overall reads

sums2 <- as.data.frame(sums1[,-2])
rownames(sums2) <- rownames(sums1)
colnames(sums2) <- "libsize"

# sums2$samples <- rownames(sums2)
# meta$samples <- meta$ID
# sums2$samples <- gsub(" ", "",sums2$samples)
# meta$samples <- gsub(" ", "",meta$samples)
# merged <- merge(sums2, meta, by = "samples")
# merged$CombinedGroup <- as.character(merged$CombinedGroup)
# lib_size_summary <- merged %>%
#   group_by(CombinedGroup) %>%  
#   summarise(
#     mean_lib_size = mean(CombinedGroup),  # replace lib_size_column with the actual column name
#     median_lib_size = median(CombinedGroup),
#     sd_lib_size = sd(CombinedGroup)
#   )


sums2 <- sums2[, ncol(sums2):1]
sums2 <- t(sums2)
View(sums2)
ncounts <- read.csv('/Users/tandre/Downloads/hires_neurons_filtered.csv')
mcounts <- read.csv('/Users/tandre/Downloads/hires_microglia_filtered.csv')
ocounts <- read.csv('/Users/tandre/Downloads/hires_oligo_filtered.csv')


#counts per million = raw counts/Library size * 1 million, so to reverse this times library size / 1 million
#multiply columns of counts by columns of libcounts 

# List of file paths and corresponding data frame names
file_paths <- c('/Users/tandre/Downloads/hires_neurons_filtered.csv',
                '/Users/tandre/Downloads/hires_microglia_filtered.csv',
                '/Users/tandre/Downloads/hires_oligo_filtered.csv')
normalized_data <- list()
for (i in 1:length(file_paths)) {

  counts <- read.csv(file_paths[i])
  norm_data <- (counts[, 3:ncol(counts)] * sums2) / 1e6
  norm_data <- round(norm_data)
  rownames(norm_data) <- counts$GeneSymbol
  normalized_data[[i]] <- norm_data
}
View(norm_data
     )
names(normalized_data) <- c('neuron_norm', 'microglia_norm', 'oligo_norm')

# neunorm <- (ncounts[,3:ncol(ncounts)]*sums2)/1e6
# rownames(neunorm) <- ncounts$GeneSymbol


#DESeq

library(DESeq2)

meta <- read.csv('/Users/tandre/Downloads/sample_sheet.csv')
meta$CombinedGroup <- paste(meta$alcohol, meta$exercise, sep = "_")
ID <- meta$ID
rownames(meta) <- gsub(" ", "", meta$ID)
meta <- meta[meta$CombinedGroup != "cage_control_cage_control", ]

# for (i in 1:length(file_paths)){
#        counts <- normalized_data[[i]]
#        colnames(counts) <- gsub(" ", "", colnames(counts))
#        normalized_data[[i]] <- counts
# }
#   

i

# neurons, should probably make this into a for loop 
View(normalized_data$neuron_norm)

meta <- meta[rownames(meta) %in% colnames(normalized_data$neuron_norm),]
normalized_data <- round(normalized_data)
View(meta)

# neurons, should probably make this into a for loop 
ddsneu <- DESeqDataSetFromMatrix(countData = normalized_data$neuron_norm,
                              colData = meta,
                              design = ~ exercise + alcohol + exercise:alcohol)
sampleDists <- dist(t(assay(ddsneu)))

sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$strain, vsd$sex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# exer x alc
dds4 <- DESeq(ddsneu,test="LRT",reduced= ~ exercise + alcohol,fitType="local")
exerc_alc <- results(dds4)
results4n <- as.data.frame(exerc_alc)[,c("stat","pvalue","padj")]
View(results4n)

#exer 

dds1 <- DESeq(ddsneu,test="LRT",reduced= ~ alcohol,fitType="local")
exerc <- results(dds1)
results1n <- as.data.frame(exerc)[,c("stat","pvalue","padj")]
View(results1n)

#alc

dds_ex_r_a <- DESeq(ddsneu,test="LRT",reduced= ~ exercise,fitType="local")
alc <- results(dds_ex_r_a)
results2n <- as.data.frame(alc)[,c("stat","pvalue","padj")]
View(results2n)
#do wald test

# microglia
ddsmicro <- DESeqDataSetFromMatrix(countData = normalized_data$microglia_norm,
                                 colData = meta,
                                 design = ~ exercise + alcohol + exercise:alcohol)
sampleDists <- dist(t(assay(ddsmicro)))

sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$strain, vsd$sex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# exer x alc
dds4 <- DESeq(ddsmicro,test="LRT",reduced= ~ exercise + alcohol,fitType="local")
exerc_alc <- results(dds4)
results4 <- as.data.frame(exerc_alc)[,c("stat","pvalue","padj")]
View(results4)

#exer 

dds1 <- DESeq(ddsmicro,test="LRT",reduced= ~ alcohol,fitType="local")
exerc <- results(dds1)
results1 <- as.data.frame(exerc)[,c("stat","pvalue","padj")]
View(results1)

#alc DEGs: 5

dds_ex_r_a <- DESeq(ddsmicro,test="LRT",reduced= ~ exercise,fitType="local")
alc <- results(dds_ex_r_a)
results2 <- as.data.frame(alc)[,c("stat","pvalue","padj")]
View(results2)

# oligodendrocytes
ddsolig <- DESeqDataSetFromMatrix(countData = normalized_data$oligo_norm ,
                                   colData = meta,
                                   design = ~ exercise + alcohol + exercise:alcohol)
sampleDists <- dist(t(assay(ddsolig)))

sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$strain, vsd$sex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# exer x alc
dds4 <- DESeq(ddsolig,test="LRT",reduced= ~ exercise + alcohol,fitType="local")
exerc_alc <- results(dds4)
results4o <- as.data.frame(exerc_alc)[,c("stat","pvalue","padj")]
View(results4o)

#exer DEGs: 3

dds1 <- DESeq(ddsolig,test="LRT",reduced= ~ alcohol,fitType="local")
exerc <- results(dds1)
results1o <- as.data.frame(exerc)[,c("stat","pvalue","padj")]
View(results1o)

#alc

dds_ex_r_a <- DESeq(ddsolig,test="LRT",reduced= ~ exercise,fitType="local")
alc <- results(dds_ex_r_a)
results2o <- as.data.frame(alc)[,c("stat","pvalue","padj")]
View(results2o)

