#Set working directory
setwd("~/")

#Read and filter methylation data for sites
methylation_data <- read.csv("Downstream_analysis/site_meth_data_hg38.txt", header = TRUE, 
                             sep = '\t', row.names = 1)
meth_sample_info <- read.csv("Annotation/annotation_file_all_samples.csv", header = TRUE, sep = ',')

methylation_data <- methylation_data[, paste('X', as.character(meth_sample_info$Sample_ID), sep = '')]
colnames(methylation_data) <- meth_sample_info$Sample_Name

#Read and filter transcriptome data for genes
transcriptome_data <- read.csv("expression_counts/normalized_gene_counts.txt", 
                               header = TRUE, sep = '\t', row.names = 1)
colnames(transcriptome_data) <- c("WT-1", "3AKO-1", "3A1-1", "3A2-1", "WT-2", "3AKO-2", "3A1-2", "3A2-2", 
                                  "WT-3", "3AKO-3", "3A1-3", "3A2-3", "WT-4", "3AKO-4", "3A1-4", "3A2-4")
common_samples <- intersect(colnames(transcriptome_data), colnames(methylation_data))
transcriptome_data <- transcriptome_data[, common_samples]
methylation_data <- methylation_data[, common_samples]

#Read list of differentially expressed genes 

rescued_genes1 <- read.table(".DESeq2_results/rescued_genesinKO_UP.txt")
rescued_genes1 <- rescued_genes1$V1
rescued_genes2 <- read.table("DESeq2_results/rescued_genesinKO_Down.txt")
rescued_genes2 <- rescued_genes2$V1
rescued_genes <- union(rescued_genes1, rescued_genes2)

#Read and filter gene and linked methylation sites info
gene_meth_sites <- read.csv("Downstream_analysis/gene_meth_sites_5000bp.txt", header = TRUE, 
                            sep = '\t')
gene_meth_sites <- subset(gene_meth_sites, gene_meth_sites$Chr != 'chrX' & gene_meth_sites$Chr != 'chrY' & 
                            gene_meth_sites$Chr != 'chrM')
rownames(gene_meth_sites) <- gene_meth_sites$Gene_id

deg_gene_list <- intersect(rescued_genes, rownames(gene_meth_sites))

deg_gene_meth_sites <- gene_meth_sites[as.character(deg_gene_list), ]
deg_gene_meth_sites <- subset(deg_gene_meth_sites, deg_gene_meth_sites$no_of_meth_sites > 0)

deg_gene_meth_sites$meth_sites <- as.character(deg_gene_meth_sites$meth_sites)

#Calculate correlation coefficient between genes and nearby methylation sites
gene_meth_site_correlation <- matrix(nrow = sum(deg_gene_meth_sites$no_of_meth_sites), ncol = 5)
rownum = 1
for (i in 1:nrow(deg_gene_meth_sites)) {
  meth_sites <- strsplit(deg_gene_meth_sites[i,8], split = ';')
  #print(meth_sites)
  for (j in 1:length(meth_sites[[1]])) {
    site_id_info <- strsplit(meth_sites[[1]][j], split = ':')
    site_id <- site_id_info[[1]][1]
    distance <- site_id_info[[1]][3]
    #print(site_id)
    corr_data_frame <- data.frame(meth=t(methylation_data[site_id,]), 
                                  expr=t(transcriptome_data[as.character(deg_gene_meth_sites[i,1]),]))
    corr_data_frame <- corr_data_frame[complete.cases(corr_data_frame),]
    rho <- cor(corr_data_frame[,1], corr_data_frame[,2], method="spearman")
    fdr <- 0
    for (k in 1:10000) {
      x <- sample(corr_data_frame[,1])
      y <- sample(corr_data_frame[,2])
      rand_rho <- cor(x, y, method="spearman")
      if(abs(rand_rho) >= abs(rho)){
        fdr <- fdr + 1
      }
    }
    fdr <- fdr/10000
    gene_meth_site_correlation[rownum, ] <- c(as.character(deg_gene_meth_sites[i,1]), site_id, distance, rho, fdr)
    rownum <- rownum + 1
  }
}
gene_meth_site_correlation <- as.data.frame(gene_meth_site_correlation)
colnames(gene_meth_site_correlation) <- c("Gene", "Site", "Distance_from_TSS", "Rho", "FDR")

gene_meth_site_correlation$Gene_name <- gene_meth_sites[as.character(gene_meth_site_correlation$Gene), "Gene_name"]
gene_meth_site_correlation <- gene_meth_site_correlation[order(gene_meth_site_correlation$FDR), ]

#Output results
write.table(gene_meth_site_correlation, file = "Downstream_analysis/rescued_genes_meth_site_correlation_5000bp_10000rep.txt", quote = FALSE,
            sep = '\t', row.names = FALSE)

#Filter for differentially methylated sites
dmp_sites <- read.csv("Differential_analysis/WT_vs_3AKO/reports/differential_methylation_data/diffMethTable_site_cmp1.csv", header = TRUE, sep = ',')
dmp_sites <- subset(dmp_sites, dmp_sites$combinedRank <= 184111)

gene_dmp_site_correlation <- subset(gene_meth_site_correlation, gene_meth_site_correlation$Site %in% dmp_sites$cgid)

write.table(gene_dmp_site_correlation, file = "Downstream_analysis/rescued_genes_dmp_site_correlation_5000bp_10000rep.txt", quote = FALSE,
            sep = '\t', row.names = FALSE)
