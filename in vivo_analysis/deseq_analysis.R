library(magrittr)
library(DESeq2)

#set working directory to project directory
setwd("~/")

#create result directory if it does not already exists
if (!dir.exists("DESeq2_results/")){
  dir.create("DESeq2_results/")
}
#read merged count data and sample info file
count_data <- read.csv("count_files/merged_gene_counts.txt", header = TRUE, sep = '\t')
colnames(count_data)[3:10] <- substr(colnames(count_data)[3:10], 1, 6)
rownames(count_data) <- count_data$Geneid
merged_counts <- count_data
count_data$Geneid <- NULL
count_data$gene_name <- NULL

col_data <- read.csv("info_files/Sequence_sample_info.csv")
count_data <- count_data[, as.character(col_data$Sample_name)]

#Make DESeq data object
dds_count <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Purity_percentage + Genotype)
dds_count <- dds_count[ rowSums(counts(dds_count)) > 1, ]
dds_count <- estimateSizeFactors(dds_count)

#Run DESeq
dds <- DESeq(dds_count, betaPrior = FALSE)
res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
res$log2FoldChange <- -1 * res$log2FoldChange
res$stat <- -1 * res$stat

#Check expression of Dnmt3a
pdf(file = "DESeq2_results/DNMT3A_gene_expression.pdf")
plotCounts(dds_count, gene = "ENSMUSG00000020661", intgroup = "Genotype")
dev.off()

#Output results to file
sink(file = "DESeq2_results/DESeq2result_summary_WT_vs_KO_controlled_by_purity.txt")
cat(summary(res))
sink()
png(file="DESeq2_results/MAplot_WT_vs_KO_controlled_by_purity.png", width = 500, height = 500)
plotMA(res, alpha=0.05, ylim=c(-5,5))
dev.off()
res_sorted <- res[order(res$padj), ]
write.table(res_sorted, file = "DESeq2_results/DESeq2result_WT_vs_KO_controlled_by_purity.txt", 
            sep="\t", quote=FALSE)

res_sorted$gene_name <- merged_counts[rownames(res_sorted), "gene_name"]
write.table(res_sorted, file = "DESeq2_results/DESeq2result_WT_vs_KO_genenames_controlled_by_purity.txt", 
            sep="\t", quote=FALSE)

