library(dplyr)

#set working directory
setwd("~/")

#Read list of in vitro rescued genes
rescued_genes1 <- read.table("invitro/DESeq2_results/rescued_genesinKO_UP.txt")
rescued_genes1 <- rescued_genes1$V1
rescued_genes2 <- read.table("invitro/DESeq2_results/rescued_genesinKO_Down.txt")
rescued_genes2 <- rescued_genes2$V1
rescued_genes <- union(rescued_genes1, rescued_genes2)

#Read list of murine genes and their overlapping methylation sites
gene_meth_sites <- read.csv("invivo/Methylation_BeadCHiP/Downstream_analysis/Baseline_analysis/gene_meth_overlap_sites_table.txt", header = TRUE, sep = '\t')

#Get orthologous human genes for murine genes
library(biomaRt)
hg_ensembl=useMart("ENSEMBL_MART_ENSEMBL")
hg_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
hg_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene"), mart = hg_mart)
hg_genes <- subset(hg_genes, hg_genes$mmusculus_homolog_ensembl_gene != '')

gene_meth_sites <- subset(gene_meth_sites, gene_meth_sites$Gene_id %in% hg_genes$mmusculus_homolog_ensembl_gene)

#Read methylation data for sites
methylation_data <- read.csv("invivo/Methylation_BeadCHiP/Downstream_analysis/Baseline_analysis/site_meth_data.txt", header = TRUE, sep = '\t')

#Read transcriptome data for genes
transcriptome_data <- read.csv("invivo/RNAseq/count_files/normalized_gene_counts.txt", header = TRUE, sep = '\t')
colnames(transcriptome_data) <- c("Dnmt3aVillin717", "Dnmt3aVillin718", "Dnmt3aVillin719", "Dnmt3aVillin720", "Dnmt3aVillin721", "Dnmt3aVillin722", "Dnmt3aVillin723", "Dnmt3aVillin724")

#Match transcriptome and methylation samples
common_samples <- intersect(colnames(transcriptome_data), colnames(methylation_data))
transcriptome_data <- transcriptome_data[, common_samples]
methylation_data <- methylation_data[, common_samples]


#Calculate correlation coefficient between genes and nearby methylation sites
gene_meth_sites$meth_id <- as.character(gene_meth_sites$meth_id)
gene_meth_sites$Gene_id <- as.character(gene_meth_sites$Gene_id)
gene_meth_site_correlation <- list()
rownum=1

for (i in 1:nrow(gene_meth_sites)) {
  corr_data_frame <- data.frame(meth=t(methylation_data[gene_meth_sites[i,]$meth_id,]), 
                                expr=t(transcriptome_data[gene_meth_sites[i,]$Gene_id,]))
  corr_data_frame <- corr_data_frame[complete.cases(corr_data_frame),]
  rho <- cor(corr_data_frame[,1], corr_data_frame[,2], method="spearman")
  fdr <- NA
  if (!is.na(rho)) {
    fdr <- 0
    for (k in 1:1000) {
      x <- sample(corr_data_frame[,1])
      y <- sample(corr_data_frame[,2])
      rand_rho <- cor(x, y, method="spearman")
      if (abs(rand_rho) >= abs(rho)) {
        fdr <- fdr + 1
      }
    }
    fdr <- fdr/1000

    
  }
  gene_meth_site_correlation[[rownum]] <- data.frame(Gene=gene_meth_sites$Gene_id[i], Site=gene_meth_sites$meth_id[i], Distance_from_TSS=gene_meth_sites$diss_from_TSS[i], 
                                                       Rho=rho, FDR=fdr)
  rownum <- rownum + 1
}
gene_meth_site_correlation <- do.call(rbind.data.frame, gene_meth_site_correlation)

gene_meth_site_correlation$Gene_name <- gene_meth_sites$Gene_name
gene_meth_site_correlation <- gene_meth_site_correlation[order(gene_meth_site_correlation$FDR), ]

#Output results
write.table(gene_meth_site_correlation, "invivo/Methylation_BeadCHiP/Downstream_analysis/Baseline_analysis/invitro_rescued_genes_site_correlation_invivo.txt", quote = FALSE, row.names = FALSE, sep = '\t')

#Compare in vivo and in vitro correlations
invitro_rescued_gene_site_correlation <- read.csv("invitro/Downstream_analysis/rescued_genes_meth_site_correlation_5000bp_10000rep.txt", sep = '\t')

gene_list <- unique(gene_meth_site_correlation[, c("Gene", "ensembl_gene_id")])
gene_list$Gene <- as.character(gene_list$Gene)
gene_list$ensembl_gene_id <- as.character(gene_list$ensembl_gene_id)
invivo_invitro_comparison <- list()

for (i in 1:nrow(gene_list)) {
  list1 <- subset(gene_meth_site_correlation, gene_meth_site_correlation$Gene == gene_list$Gene[i])
  list1 <- subset(list1, !is.na(list1$Rho))
  list2 <- subset(invitro_rescued_gene_site_correlation, invitro_rescued_gene_site_correlation$Gene == gene_list$ensembl_gene_id[i])
  tup <- c(gene_list$Gene[i], gene_list$ensembl_gene_id[i], list1[which(abs(list1$Rho) == max(abs(list1$Rho)))[1], "Rho"], list1[which(abs(list1$Rho) == max(abs(list1$Rho)))[1], "FDR"], 
           list2[which(abs(list2$Rho) == max(abs(list2$Rho)))[1], "Rho"], list2[which(abs(list2$Rho) == max(abs(list2$Rho)))[1], "FDR"])
  invivo_invitro_comparison[[i]] <- tup
  
}

invivo_invitro_comparison <- do.call(rbind.data.frame, invivo_invitro_comparison)
colnames(invivo_invitro_comparison) <- c("Gene.invivo", "Gene.invitro", "Rho.invivo", "FDR.invivo", "Rho.invitro", "FDR.invitro")
invivo_invitro_comparison$Rho.invivo <- as.numeric(as.character(invivo_invitro_comparison$Rho.invivo))
invivo_invitro_comparison$Rho.invitro <- as.numeric(as.character(invivo_invitro_comparison$Rho.invitro))
invivo_invitro_comparison$FDR.invivo <- as.numeric(as.character(invivo_invitro_comparison$FDR.invivo))
invivo_invitro_comparison$FDR.invitro <- as.numeric(as.character(invivo_invitro_comparison$FDR.invitro))

write.table(invivo_invitro_comparison, "invivo/Methylation_BeadCHiP/Downstream_analysis/invitro_rescued_genes_correlation_invivo_invitro_comparison.txt", quote = FALSE, row.names = FALSE, sep = '\t')

#Plot in vivo-in vitro correlation comparison
invivo_invitro_comparison <- subset(invivo_invitro_comparison, invivo_invitro_comparison$FDR.invivo < 0.05 | invivo_invitro_comparison$FDR.invitro < 0.05)

invivo_invitro_comparison$sig <- ifelse(invivo_invitro_comparison$FDR.invivo < 0.05 & invivo_invitro_comparison$FDR.invitro < 0.05, "both",
                                        ifelse(invivo_invitro_comparison$FDR.invitro < 0.05, "invitro", "invivo"))
invivo_invitro_comparison[is.na(invivo_invitro_comparison$sig), "sig"] <- "invitro"

line_colors <- c("#714C02", "#01587A", "#024E37")
fill_colors <- c("#9D6C06", "#077DAA", "#026D4E")

pdf("invivo/Methylation_BeadCHiP//Downstream_analysis/invitro_rescued_genes_correlation_invivo_invitro_comparison.pdf")
p <- ggplot()
p <- p + geom_point(data=subset(invivo_invitro_comparison, invivo_invitro_comparison$sig=="invitro"), mapping = aes(x=Rho.invitro, y=Rho.invivo, color="#01587A", fill="#077DAA"), 
                    size=3, shape=21, alpha=0.5)
p <- p + geom_point(data=subset(invivo_invitro_comparison, invivo_invitro_comparison$sig=="invivo"), mapping = aes(x=Rho.invitro, y=Rho.invivo, color="#024E37", fill="#026D4E"), 
                    size=3, shape=21, alpha=0.5)
p <- p + geom_point(data=subset(invivo_invitro_comparison, invivo_invitro_comparison$sig=="both"), mapping = aes(x=Rho.invitro, y=Rho.invivo, color="#714C02", fill="#9D6C06"), 
                    size=3, shape=21, alpha=0.7)
p <- p + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p <- p + scale_fill_identity(guide='legend', labels=c("invivo", "invitro", "both")) + 
  scale_color_manual(values = c("#01587A"="#01587A", "#024E37"="#024E37","#714C02"="#714C02"))
p <- p + xlab("invitro correlation") + ylab("invivo correlation")
p <- p + theme_minimal() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=14), 
                                 legend.title = element_text(size=0), plot.title = element_text(size=16, face="bold", hjust = 0.5), 
                                 axis.line = element_line(size=0.7, colour = "black")) 
p
dev.off()
