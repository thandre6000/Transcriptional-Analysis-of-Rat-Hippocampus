
library(tidyverse)
library(Seurat)
library(dplyr)

#Seurat object name from this set is "ss.seurat"

load("/Users/theoandre/Downloads/Seurat.ss.rda")

#filtered to include astro, micro pvm, and the like because (I believe) the HIP 
#region label only includes specialized hippocampal cells. More info is available at 
#https://star-protocols.cell.com/protocols/1386#sec1.2

sc_data<-subset(x = ss.seurat,
                subset = region_label == "HIP" | subclass_label %in% c("Astro","Micro-PVM","Endo","Oligo","VLMC"))
sc_data<-subset(x= sc_data,
                subset = subclass_label %in% c("Astro","CA1-ProS","CA3","DG","Endo","Lamp5","Oligo","Micro-PVM","Pvalb","Sncg","Sst","Vip"))
#set as factor
sc_data$subclass_label <- factor(sc_data$subclass_label)
#drop levels
sc_data$subclass_label<- droplevels(x = sc_data$subclass_label)
#rename cell types
Idents(sc_data) <- "subclass_label"
sc_data <- RenameIdents(sc_data,
                     "Astro" = "Astro",
                     "CA1-ProS" = "Pyr",
                     "CA3" = "Pyr",
                     "Endo" = "Endo",
                     "DG" = "Dent",
                     "Lamp5" = "Inter",
                     "Micro-PVM" = "Micro",
                     "Oligo" = "Oligo",
                     "Pvalb" = "Inter",
                     "Sncg" = "Inter",
                     "Sst" = "Inter",
                     "Vip" = "Inter")

Idents(sc_data)<-factor(Idents(sc_data), levels = sort(levels(sc_data)))

#CPM normalization
sc_data <- NormalizeData(sc_data, normalization.method = "RC", scale.factor = 1e6)

#rename all individual cells to their respective cell type

#append subclass as the first row of the df

subclass <- Idents(sc_data)
counts <- as.data.frame(sc_data@assays$RNA@counts)
subclass <- as.data.frame(subclass)
subclass <- t(subclass)
countsnamed <- rbind(subclass, counts)
subclassfact <- as.numeric(as.factor(subclass))
countsnamed[1,]<- subclassfact
colnames(countsnamed)<- NULL

write.table(countsnamed, file = "/Users/theoandre/Desktop/Lab/Deconvolution/allencounts.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#test to see
write.table(head(countsnamed), file = "/Users/theoandre/Desktop/Lab/Deconvolution/reference2.tsv", row.names = TRUE)







#if you want to explore a seurat object, you can use assays$RNA@
#export the tsv- done, rest is just for reference
write.table(sc_data@assays$RNA@counts, file = "/Users/theoandre/Desktop/DeconvolutionReferencetable1.tsv", row.names = TRUE)

#export as tsv/normalize counts (can also remove low read counts)
install.packages("scTenifoldNet")
library(scTenifoldNet)

ratsOutput <- cpmNormalization(rats)
View(ratsOutput)

# Visualizing the differences
oldPar <- par(no.readonly = TRUE)

par(
  mfrow = c(1, 2),
  mar = c(3, 3, 1, 1),
  mgp = c(1.5, 0.5, 0)
)
plot(
  Matrix::colSums(rats),
  ylab = 'Library Size',
  xlab = 'Cell',
  main = 'Before CPM Normalization'
)
plot(
  Matrix::colSums(ratsOutput),
  ylab = 'Library Size',
  xlab = 'Cell',
  main = 'After CPM Normalization'
)

par(oldPar)

