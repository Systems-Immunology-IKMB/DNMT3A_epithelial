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
get_group_go_data <- function(groups, n_terms){
  plist <- vector('list', length(groups))
  go_list <- c()
  go_terms <- c()
  for (i in 1:length(groups)) {
    #print(modules[i])
    go_results <- read.csv(groups[i], 
                           header = TRUE, sep = '\t')
    go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
    go_results <- subset(go_results, go_results$Significant > 1)
    #print(nrow(go_results))
    top_go_results <- filter_go_results(go_results, n_terms)
    #print(top_go_results)
    go_list <- c(go_list, as.character(top_go_results$GO.ID))
    go_terms <- c(go_terms, as.character(top_go_results$Term2))
    plist[[i]] <- go_results
    
  }
  
  p_vector <- c()
  n_vector <- c()
  for (i in 1:length(groups)) {
    p_vector <- c(p_vector, paste("p", i, sep = '_'))
    n_vector <- c(n_vector, paste("n", i, sep = '_'))
  }
  
  go_data <- data.frame(GO.ID=go_list, Term=go_terms)
  go_data[, p_vector] <- NA
  go_data[, n_vector] <- NA
  
  for (i in 1:length(groups)) {
    df <- plist[[i]]
    rownames(df) <- as.character(df$GO.ID)
    group_go_data <- df[as.character(go_data$GO.ID),]
    go_data[, (i+2)] <- group_go_data$Fisher.elim
    go_data[, (i+2+length(groups))] <- group_go_data$Significant/group_go_data$Annotated
  }
  
  return(go_data)
}


groups <- c("invitro/AF_Caco2_DNMT3A/DESeq2_results/WT_vs_3AKO/DESeq2result_WT_vs_3AKO_unpaired_GO_up.txt", 
            "invivo/RNAseq/DESeq2_results/DESeq2result_WT_vs_KO_up_GO.txt")
group_go_data_up <- get_group_go_data(groups, 20)

#Reformat data for plotting
p_vector <- colnames(group_go_data_up)[grep("p_", colnames(group_go_data_up))]
n_vector <- colnames(group_go_data_up)[grep("n_", colnames(group_go_data_up))]
plot_data_up <- melt(group_go_data_up[, c("GO.ID", "Term", p_vector)], id.vars = c("GO.ID","Term"))
plot_data_up$n <- unlist(group_go_data_up[, n_vector])
plot_data_up$Group <- groups[as.numeric(sub("p_", "", plot_data_up$variable))]
plot_data_up$value <- as.numeric(as.character(plot_data_up$value))
plot_data_up$p <- -1*log10(plot_data_up$value)
plot_data_up$Term2 <- reorder(plot_data_up$Term, 1:nrow(plot_data_up))

p_vector <- colnames(group_go_data_down)[grep("p_", colnames(group_go_data_down))]
n_vector <- colnames(group_go_data_down)[grep("n_", colnames(group_go_data_down))]
plot_data_down <- melt(group_go_data_down[, c("GO.ID", "Term", p_vector)], id.vars = c("GO.ID","Term"))
plot_data_down$n <- unlist(group_go_data_down[, n_vector])
plot_data_down$Group <- groups[as.numeric(sub("p_", "", plot_data_down$variable))]
plot_data_down$value <- as.numeric(as.character(plot_data_down$value))
plot_data_down$p <- -1*log10(plot_data_down$value)
plot_data_down$Term2 <- reorder(plot_data_down$Term, 1:nrow(plot_data_down))
plot_data_down$p <- -1 * plot_data_down$p

plot_data <- rbind(plot_data_up, plot_data_down)
plot_data <- plot_data_up
plot_data$Term2 <- str_wrap(plot_data$Term2, width = 60)
plot_data$Term2 <- reorder(plot_data$Term2, 1:nrow(plot_data))

#Plot result
cairo_ps("invivo/RNAseq/DESeq2_results/invitro_invivo_top_GO_terms.ps", width = 9, height = 14)
p <- ggplot(plot_data, mapping = aes(x=Group, y=Term2, size=n, color=p))
p <- p + geom_point()
p <- p + scale_x_discrete(limits=groups)
p <- p + xlab("") + ylab("GO Term")
p <- p + scale_colour_gradient2(high="#660000", low="#000066", mid="white", name = "-log10p") 
p <- p + scale_size(name = "Gene-ratio")
p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                             axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1), 
                             axis.title = element_text(size = 20))
p
dev.off()


