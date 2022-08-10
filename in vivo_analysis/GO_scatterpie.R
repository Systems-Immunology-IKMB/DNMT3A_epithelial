library(GO.db)
library(stringr)
library(ggplot2)
library(reshape2)
#Set working directory
setwd("~/")

#Filter GO results for significance, minimum number of genes, ontology, unique gene sets and maximum number of terms
filter_go_results <- function(go_results, top_num){
  go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
  go_results <- subset(go_results, go_results$Significant > 1)
  go_results <- subset(go_results, go_results$ont == "BP")
  if (nrow(go_results) > 0) {
    go_results <- go_results[order(go_results$Fisher.elim), ]
    #print(nrow(go_results))
    exclude_indices <- c()
    unique_gene_sets <- c()
    for(i in 1:length(go_results$GO.ID)){
      #print(go_results_sorted$Genes[i])
      if(as.vector(go_results$Genes[i]) %in% unique_gene_sets){
        exclude_indices <- c(exclude_indices, i)
      }else{
        unique_gene_sets <- c(unique_gene_sets, as.vector(go_results$Genes[i]))
      }
    }
    if(length(exclude_indices) > 0){
      go_results_filtered <- go_results[-exclude_indices, ]
    }else{
      go_results_filtered <- go_results
    }
    #go_results_filtered <- go_results
    if(length(go_results_filtered$GO.ID) > top_num){
      top <- go_results_filtered[1:top_num,]
    }else{
      top <- go_results_filtered
    }
    
    top$p <- -1*log10(top$Fisher.elim)
    top$Term <- str_wrap(Term(as.vector(top$GO.ID)), width = 40)
    top$Term <- Term(as.vector(top$GO.ID))
    top$Term2 <- reorder(top$Term, top$p)
    
    return(top)
  }
  
}
#Group GO results 
get_group_go_data <- function(groups, n_terms, dmp_data){
  plist <- vector('list', length(groups))
  go_list <- c()
  go_terms <- c()
  for (i in 1:length(groups)) {
   
    go_results <- read.csv(groups[i], 
                           header = TRUE, sep = '\t')
    go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
    go_results <- subset(go_results, go_results$Significant > 1)
    go_results <- get_hyper_hypo_number(go_results, dmp_data[[i]])
   
    top_go_results <- filter_go_results(go_results, n_terms)
   
    go_list <- c(go_list, as.character(top_go_results$GO.ID))
    go_terms <- c(go_terms, as.character(top_go_results$Term2))
    plist[[i]] <- go_results
    
  }
 
  go_data <- data.frame()
  
  for (i in 1:length(groups)) {
   
    df <- plist[[i]]
    rownames(df) <- as.character(df$GO.ID)
   
    group_go_data <- df[intersect(go_list, rownames(df)),]
    group_go_data$group <- groups[i]
   
    go_data <- rbind(go_data, group_go_data)
    
  }
  go_data$n <- go_data$Significant/go_data$Annotated
  return(go_data)
}

#Get number of hyper and hypomethylated genes for each GO term
get_hyper_hypo_number <- function(go_results, dmpr_results){
  hypermethylated_promoters <- subset(dmpr_results, dmpr_results$mean.mean.diff < 0)$symbol
  hypomethylated_promoters <- subset(dmpr_results, dmpr_results$mean.mean.diff > 0)$symbol
  go_results$Genes <- gsub("^,", "", go_results$Genes)
  go_results$Genes_vector <- sapply(go_results$Genes, strsplit, ',')
  
  f <- function(x, y){return(length(intersect(x,y)))}
  
  go_results$Hypermethylated <- sapply(go_results$Genes_vector, f, hypermethylated_promoters)
  go_results$Hypomethylated <- sapply(go_results$Genes_vector, f, hypomethylated_promoters)
  
  return(go_results)
}

groups <- c("Downstream_analysis/WT_DSS5_analysis/DMPr_DSS5_WT_only_GO.txt", 
            "Downstream_analysis/KO_DSS5_analysis/DMPr_DSS5_KO_only_GO.txt", 
            "Downstream_analysis/WT_DSS5_analysis/DMPr_DSS5_WT_KO_overlap_GO.txt", 
            "Downstream_analysis/WT_DSS12_analysis/DMPr_DSS12_WT_only_GO.txt", 
            "Downstream_analysis/KO_DSS12_analysis/DMPr_DSS12_KO_only_GO.txt", 
            "Downstream_analysis/WT_DSS12_analysis/DMPr_DSS12_WT_KO_overlap_GO.txt")
dmpr_results <- list(read.csv("Analysis/reports_WT_DSS5/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"), 
                     read.csv("Analysis/reports_KO_DSS5/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"), 
                     read.csv("Analysis/reports_WT_DSS5/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"), 
                     read.csv("Analysis/reports_WT_DSS12/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"), 
                     read.csv("Analysis/reports_KO_DSS12/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"), 
                     read.csv("Analysis/reports_WT_DSS12/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"))
group_go_data <- get_group_go_data(groups, 50, dmpr_results)

group_go_data$group <- gsub("Downstream_analysis/WT_DSS5_analysis/DMPr_", "", group_go_data$group)
group_go_data$group <- gsub("Downstream_analysis/WT_DSS12_analysis/DMPr_", "", group_go_data$group)
group_go_data$group <- gsub("Downstream_analysis/KO_DSS5_analysis/DMPr_", "", group_go_data$group)
group_go_data$group <- gsub("Downstream_analysis/KO_DSS12_analysis/DMPr_", "", group_go_data$group)
group_go_data$group <- gsub("_GO.txt", "", group_go_data$group)
group_go_data$group <- gsub("_", " ", group_go_data$group)

plot_data <- group_go_data
plot_data$p <- -1*log10(plot_data$Fisher.elim)
plot_data$Term <- Term(as.vector(plot_data$GO.ID))
plot_data$Term2 <- str_wrap(plot_data$Term, width = 60)
plot_data$Term2 <- reorder(plot_data$Term2, 1:nrow(plot_data))

plot_data <- transform(plot_data, y=match(GO.ID, unique(GO.ID)))
plot_data <- transform(plot_data, x=match(group, unique(group)))
plot_data$radius <- plot_data$p * 0.2

#Plot results
pdf("Downstream_analysis/DSS_promoter_selected_GO_scatterpie.pdf", width = 8, height = 9)
p <- ggplot() + geom_scatterpie(aes(x=x, y=y, r=radius), data = plot_data, cols = c("Hypermethylated", "Hypomethylated"), color=NA) + coord_equal()
p <- p + geom_point()
p <- p + scale_x_continuous(breaks = 1:6, labels = unique(selected_plot_data$group))
p <- p + scale_y_continuous(breaks = 1:max(selected_plot_data$y), labels = unique(selected_plot_data$Term2))
p <- p + xlab("") + ylab("GO Term")
p <- p + scale_fill_manual(values = c('#990000', '#004C99'))
p <- p + geom_scatterpie_legend(selected_plot_data$radius, x = 1, y = 23)
p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                            axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1), 
                            axis.title = element_text(size = 20))
p
dev.off()
