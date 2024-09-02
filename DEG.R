getwd()
setwd('C:\\Users\\chaud\\Pictures\\edgene')

count <- as.matrix(read.csv("count_data.csv", sep = ",", row.names = "gene_id"))
head(count)

meta <- read.csv("metadata.csv", row.names = 1)
head(meta)

all(rownames(meta) == colnames(count))

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~condition)
dim(dds)

keep <- rowSums(counts(dds))>=10  #remove low count genes
dds <- dds[keep,]
dim(dds)

#normalisation using vst
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
head(assay(vsd),3)

res <- results(dds, contrast = c("condition", "treated", "untreated"), alpha=0.05)
summary(res)
write.csv(res, "DEG.csv")
