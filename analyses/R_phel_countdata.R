#Archivo para generar datos de expresion con el archivo
#Phel_countdata.txt"
## try http if https is not available

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library(DESeq2)
../../mdelrio2-btae-nb/analyses/Phel_countdata.txt
#C:/users/curso06/Desktop/mdelrio2-btae-nb/Phel_countdata.txt"
#directorio por omision "Libraries\Documents" en esta maquina, es necesario
#verificar en la maquina con la que se este trabajando.

data <- read.table("c:Phel_countdata.txt", header = T, sep = "\t")
rownames(data) <- data$Feature
data <- data[,-1]
# Build objets
# Specify which columns are in which groups
deseq2.colData <- data.frame(condition=factor(c(rep("Treated", 3), rep("Control", 3))), 
                             type=factor(rep("single-read", 6)))
rownames(deseq2.colData) <- colnames(data)
deseq2.dds <- DESeqDataSetFromMatrix(countData = data,
                                     colData = deseq2.colData, 
                                     design = ~ condition)


# Run Analysis
deseq2.dds <- DESeq(deseq2.dds)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]
head(deseq2.res)
# Count number of hits with adjusted p-value less then 0.05
dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])


tmp <- deseq2.res
# The main plot
plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="darkgray",
     main="DEG Virus Exposure  (pval <= 0.05)",
     xlab="mean of normalized counts",
     ylab="Log2 Fold Change")
# Getting the significant points and plotting them again so they're a different color
tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="red")
# 2 FC lines
abline(h=c(-1,1), col="blue")
write.table(tmp.sig, "c:Phel_DEGlist.tab", row.names = T)

head("c:Phel_DEGlist.tab")
