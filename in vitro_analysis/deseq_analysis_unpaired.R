library(magrittr)
library(DESeq2)

#set working directory to project directory
setwd("/home/user/Projects/AF_Caco2_DNMT3A/")
#create result directory if it does not already exists
if (!dir.exists("DESeq2_results/")){
  dir.create("DESeq2_results/")
}

species <- "human"

#define conditions to compare
condition1 <- 'siWT'
condition2 <- '3AKO'

#read merged count data and sample info file
count_data1 <- read.table("expression_counts/merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
count_data2 <- read.table("expression_counts/merged_gene_counts_si_control.txt", header = TRUE, sep = '\t', row.names = 1)
count_data <- cbind(count_data1, count_data2)
col_data <- read.csv("info_files/Sample_info.csv")
colnames(count_data) <- col_data$Sample.name

#filter col_data and count_data based on conditions
col_data <- subset(col_data, col_data$Condition %in% c(condition1, condition2))
#col_data$Condition <- as.character(col_data$Condition)
#col_data$Sample_id <- as.character(col_data$Sample_id)
#col_data$Sample.name <- as.character(col_data$Sample.name)
#col_data$Barcode <- as.character(col_data$Barcode)
count_data <- count_data[, as.character(col_data$Sample.name)]

#create output directory if it does not already exists
output_directory <- paste("DESeq2_results/", condition1, "_vs_", sub('%', '', sub('\\.', '', condition2)), "/", sep = '')
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

#Make DESeq data object
dds_count <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data[,c("Sample.name", "Condition")], design = ~ Condition)
dds_count <- dds_count[ rowSums(counts(dds_count)) > 1, ]
dds_count <- estimateSizeFactors(dds_count)

#Run DESeq
dds <- DESeq(dds_count, betaPrior = FALSE)
res <- results(dds, independentFiltering = TRUE, alpha = 0.05)

#Output results to file
sink(file = paste(output_directory, "DESeq2result_summary_", condition1, "_vs_", sub('%', '', sub('\\.', '', condition2)), "_unpaired.txt", sep = ""))
cat(summary(res))
sink()
png(file=paste(output_directory, "MAplot_", condition1, "_vs_", sub('%', '', sub('\\.', '', condition2)), "_unpaired.png", sep = ""), width = 500, height = 500)
plotMA(res, alpha=0.05, main=paste(condition1, "vs.", condition2, sep = " "), ylim=c(-7,7))
dev.off()
res_sorted <- res[order(res$padj), ]
write.table(res_sorted, file = paste(output_directory, "DESeq2result_", condition1, "_vs_", sub('%', '', sub('\\.', '', condition2)), "_unpaired.txt", sep = ""), 
            sep="\t", quote=FALSE)

#Add gene symbols
if (species == "mouse") {
  mart_export <- read.table("/home/user/Projects/reference/mouse/GRCm38_mart_export.txt", header = T, sep="\t")
}else if (species == "human"){
  mart_export <- read.table("/home/user/Projects/reference/human/GRCh38_mart_export.txt", header = T, sep="\t") 
}

#mart_export <- read.table("mart_export.txt", header = T, sep = '\t')
unique_mart <- subset(mart_export, duplicated(mart_export$Gene_stable_ID) == FALSE)
rownames(unique_mart) <- unique_mart$Gene_stable_ID
res_sorted$gene <- unique_mart[rownames(res_sorted), ]$Gene_name
write.table(res_sorted, file = paste(output_directory, "DESeq2result_genenames_", condition1, "vs", sub('%', '', sub('\\.', '', condition2)), "_unpaired.txt", sep = ""), 
            sep="\t", quote=FALSE)
