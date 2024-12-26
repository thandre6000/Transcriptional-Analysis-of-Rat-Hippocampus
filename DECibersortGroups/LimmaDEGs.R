



BiocManager::install("limma")

library(limma)

logCM <- read.csv('/Users/tandre/Desktop/hires_microglia_filtered_logcpm.csv')
rownames(logCM) <- logCM$GeneSymbol
logCM <- logCM[,-1]
#format for limma
logCMt <- t(logCM)



meta <- read.csv('/Users/tandre/Downloads/sample_sheet.csv')
meta$CombinedGroup <- paste(meta$alcohol, meta$exercise, sep = "_")
ID <- meta$ID
meta <- as.data.frame(meta[,c(5),drop = FALSE])
rownames(meta) <- ID

View(design)

design <- model.matrix(~ 0 + CombinedGroup, data = meta)
colnames(design) <- c("binge_ex", "binge_sed", "cc", "control_ex", "control_sed")
design <- design[,-3]
# design <- as.matrix(meta)
rownames(design) <- gsub(" ", "", rownames(design))

design1 <- design[rownames(logCMt),]
fit <- lmFit(logCM, design1)
# fit<- eBayes(fit, trend=TRUE)
contrasts <- makeContrasts(cond1 = binge_ex-control_sed, cond2 = binge_ex-control_ex, levels=design1)

fit_contrasts <- contrasts.fit(fit, contrasts)
fit_contrasts <- eBayes(fit_contrasts)

genes <-topTable(fit_contrasts, coef= "cond1", number= nrow(logCM))
View(genes)
genesnames<-rownames(genes)
fit_contrasts <- contrasts.fit(fit, contrasts)


#heatmap

install.packages("pheatmap"
                 )
library(pheatmap)
group_means <- t(logCMt) %*% design1/ colSums(design)
top_genes <- apply(group_means, 1, var)
selected_genes <- names(sort(top_genes, decreasing = TRUE))[1:50]
group_means_selected <- group_means[selected_genes, ]

annotation_colors <- list(
  Group = c("binge_ex" = "lightblue", 
            "binge_sed" = "lightgreen", 
            "control_ex" = "lightpink", 
            "control_sed" = "salmon")
)

annotation <- data.frame(Group = factor(colnames(group_means_selected)))
# Create the heatmap
png("/Users/tandre/Desktop/heatmap.png", width = 800, height = 600)
pheatmap(group_means,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100))
