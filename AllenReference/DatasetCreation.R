---
  title: "Reference Data Set Creation"
author: "Theo Andre"
date: "2024-03-31"
---
  
  #set wd and download packages
  

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
install.packages('data.table')

# load necessary packages
library(tidyverse)
library(data.table)

#Read in reference data

counts <- fread("/project/leasure/singlecell/matrix.csv", data.table = FALSE)
metadata <- read.csv(file = "/Users/theoandre/Downloads/metadata.csv", row.names = 1)

rownames(counts) <- counts$sample_name
counts <- counts[,-!names(counts) %in% c("sample_name")]
counts <- as.matrix(counts)

transposed_counts <- t(counts)
rm(counts)
gc()

#check that each column is backcode and each row is gene (will take a while)

View(transposed_counts)

#build the object

library(Seurat)
sc_data<- CreateSeuratObject(counts = transposed_counts,
                             meta.data = metadata, min.cells = 0, min.features = 0,
                             project ="Alcoholic Rats")

rm(transposed_counts)
rm(metadata)
gc()

sc_data<-subset(x = sc_data,
                subset = region_label == "HIP")
sc_data<-subset(x= sc_data,
                subset = subclass_label=="Astro","CA1-ProS","CA3","DG","Endo","Lamp5","Oligo","Micro-PVM","Pvalb","Sncg","Sst","Vip")

sc_data$subclass_label<- droplevels(x = sc_data$subclass_label)
Idents(sc_data) <- "subclass_label"
c_data<-RenameIdents(sc_data,
                     "Astro" = "Astrocytes",
                     "CA1-ProS" = "Pyramidal",
                     "CA3" = "Pyramidal",
                     "Endo" = "Endothelial",
                     "DG" = "Dentate"
                     "Lamp5" = "Interneurons",
                     "Micro-PVM" = "Microglia",
                     "Oligo" = "Oligodendrocytes",
                     "Pvalb" = "Interneurons",
                     "Sncg" = "Interneurons",
                     "Sst" = "Interneurons",
                     "Vip" = "Interneurons",
                     )
Idents(sc_data)<-factor(Idents(sc_data), levels = sort(levels(sc_data)))
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
sc_data <- ScaleData(sc_data)

#visualize

sc_data <- RunPCA(sc_data)
ElbowPlot(sc_data, ndims = 50)

sc_data <- RunTSNE(sc_data, dims = 1:20)
sc_data <- RunUMAP(sc_data, dims = 1:20)

library(cowplot)
library(ggplot2)
library(patchwork)

cols <-c("limegreen", #Astrocytes
         
         "steelblue", #Endothelial
         
         "mediumorchid4", #Interneurons
         
         "firebricks2", #Microglia
         
         "magenta", #Mural
         
         "gray52", #Oligodendrocytes
         
         "tan1") #Pyramidal

pca_plot<-DimPlot(sc_data, reduction = "pca", pt.size = 0.1, label = T, cols = cols)
umap_plot<-DimPlot(sc_data, reduction = "umap", pt.size = 0.1, label = T, cols = cols)

legend<-get_legend(umap_plot)
layout<-c (area(1, 1, 1, 2),#PCA
           
           area(1, 5, 1, 6),#UMAP
           
           area(1, 7))#Legend

fig1 <- pca_plot + tsne_plot + umap_plot + legend + plot_layout(design = layout) & NoLegend()
ggsave(filename = "Fig_1.png", plot = fig1, width = 13, height = 3.75, dpi = 300)


#write sc data to a tsv, keeping only layers

write.table(sc_data@assays$RNA@layers, file = "/Users/theoandre/Desktop/Deconvolution/Referencetable.tsv", row.names = TRUE)


#this method did not work
#Counts <- GetAssayData(sc_data)
#write.table(Counts, file = "/Users/theoandre/Desktop/Deconvolution/Referencetable.tsv", row.names = TRUE)

#you can veiw this in R
#save

saveRDS(sc_data, "Single-cell_custom_subset.rds")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")

library(rhdf5)
h5ls("expression_matrix.hdf5")
matrix <- h5read("/Users/theoandre/Downloads/expression_matrix.hdf5")




# Calculate the number of cells with 0
zero_cells <- sum(counts_numeric_part == 0, na.rm = TRUE)

# Calculate total number of cells in the numeric part
total_cells_numeric_part <- length(counts_numeric_part) * nrow(counts_numeric_part)

# Calculate percentage of cells with 0
percentage_zero <- (zero_cells / total_cells_numeric_part) * 100
percentage_zero

