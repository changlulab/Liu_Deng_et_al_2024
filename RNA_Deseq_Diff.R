M3_RNA_count<-read.csv("E:/R files/M3_RNA_count.csv",row.names = 1)
M3_RNA_count<-as.matrix(M3_RNA_count)
M3_RNA_count_coldata<-read.csv("E:/R files/M3_RNA_count_coldata.csv",row.names = 1)
M3_RNA_count_coldata<-as.matrix(M3_RNA_count_coldata)
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = M3_RNA_count, colData =M3_RNA_count_coldata, design = ~ condition)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, file = "E:/R files/M3_RNA_count_Diff.csv")