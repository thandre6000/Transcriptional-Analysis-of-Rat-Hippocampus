samples <- read.csv("/Users/theoandre/Downloads/sample_sheet.csv")
samples$condition <- paste(samples$exercise, samples$alcohol)
samples <- select(samples, c("ID", "condition"))
samples$condition <- replace(samples$condition, samples$condition == "cage_control cage_control", "cage_control")

colnames(samples)[1] <- "gene_id"
colnames <- samples$gene_id

samplestrans <- t(samples)
colnames(samplestrans)<- colnames
samplestrans <- samplestrans[!rownames(samplestrans) == "gene_id", , drop = FALSE]