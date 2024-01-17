## Script for doing differential expression analysis using DESeq2 package in R

rm(list=ls())

#install packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("ggplot2")

#upload required libraries

library(DESeq2)
library(ggplot2)

## Data processing
data <- read.table("Counts.txt", sep="\t", header=TRUE, check.names=FALSE)
row.names(data) <- make.names(data$Gene, unique = TRUE, allow_ = FALSE)

## DESeq2 ##
## Design matrix
design <- read.table(
        file="designMatrix.txt",
        header=TRUE,
        sep="\t",
        row.names=1
        )

## Differential expression analysis (Untreated vs Treated)
this_analysis <- "Diseased_vs_Control"
this_data <- data[,c(2:6)]
this_design <- design[c(1:5),]
dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ Condition
		)
dds <- dds[rowSums(counts(dds)) > 5,]	## Remove the low read count genes
dds <- DESeq(dds)
res <- results(dds,contrast=c("Condition", "Diseased", "Control"))
	
pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
plotDispEsts( dds )
plotMA(res)
dev.off()
	
## Writing outputs		
res <- as.data.frame(res)
res$Gene <- row.names(res)	
res$padj.BH <- p.adjust(res$pvalue, method="BH") 
res <- merge(res, data, by="Gene")
res <- merge(res, tpm, by="Gene")
res <- res[order(res$padj.BH),]
write.table(res, file = paste(this_analysis, "_deseq2_results_all_genes.txt", sep=""), sep = "\t", row.names=F, quote=F)
res_sig <- subset(res, res$padj.BH < 0.05 & (res$log2FoldChange > 0.5 | res$log2FoldChange < -0.5))
write.table(res_sig, file = paste(this_analysis, "_deseq2_results_significant_genes.txt", sep=""), sep = "\t", row.names=F, quote=F)
write.table(res_sig[,c(1,3)], file = paste(this_analysis, "_deseq2_results_significant_genes.rnk", sep=""), sep = "\t", row.names=F, quote=F, col.names=F)
write.table(res_sig[,c(1,17:22)], file = paste(this_analysis, "_input_for_heatmap.txt", sep=""), sep = "\t", row.names=F, quote=F, col.names=T)

## Volcano plot ##
pdf(paste(this_analysis, "_volcano.pdf", sep=""), width=8, height=8)

with(data, plot(logFC, -log10(adj.P.Val), pch=20, ylim=c(1,3), cex=0.5, col="red", xlab="∇β-value", ylab="-Log10(Adj.Pvalue) ", cex.lab=1.5, cex.axis=1.5, cex.names=1.5))
with(subset(data, adj.P.Val > 0.01 | (data$logFC < 0)), points(logFC, -log10(adj.P.Val), pch=20, col="green", cex=0.5))
with(subset(data, adj.P.Val > 0.01 | (data$logFC < 0.2 & data$logFC > -0.2)), points(logFC, -log10(adj.P.Val), pch=20, col="grey", cex=0.5))
legend("topright", legend=c("Downregulated DEGs", "Upregulated DEGs"), pch=20, col=c("green", "red"), cex=0.7)
abline(h=2, col = "blue",lty=2)
abline(v=0.2, col = "blue",lty=2)
abline(v=-0.2, col = "blue",lty=2)
dev.off()
