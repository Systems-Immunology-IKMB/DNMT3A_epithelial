library(topGO)
library(DESeq2)
library(genefilter)
library(geneplotter)
library(plyr)
library(GO.db)
library(org.Mm.eg.db)

#Set working directory
setwd("~/")

#Read gene lists
dds_res <- read.table("DESeq2_results/DESeq2result_WT_vs_KO_genenames.txt", 
                      header = TRUE, sep = '\t')
degs_tab <- subset(dds_res, dds_res$padj < 0.05)
up_degs <- rownames(subset(degs_tab, degs_tab$log2FoldChange > 0))
down_degs <- rownames(subset(degs_tab, degs_tab$log2FoldChange < 0))

overallBaseMean <- as.matrix(dds_res[, "baseMean", drop = F])

#Identify genes with similar gene expression level
sig_idx <- match(up_degs, rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, up_degs)

multidensity( list( 
  all= log2(dds_res[,"baseMean"]+1) ,
  foreground =log2(dds_res[up_degs, "baseMean"]), 
  background =log2(dds_res[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

#GO enrichment analysis
onts = c( "MF", "BP", "CC" )

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(up_degs,  backG) 
inSelection =  geneIDs %in% up_degs 
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]

tab = as.list(onts)
names(tab) = onts
for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  ## look at results
  if(length(nodes(graph(tgd))) < 200){
    tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                          Fisher.classic = resultTopGO.classic,
                          orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
  }else{
    tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.elim" , topNodes = 200)}
  
}

tab$MF$ont <- "MF"
tab$BP$ont <- "BP"
tab$CC$ont <- "CC"
topGOResults <- rbind.fill(tab)
go_results <- topGOResults

#Add gene names
for(i in 1:length(go_results$GO.ID)){
  go_id <- as.vector(go_results[i,1])
  alleges <- get(go_id, org.Mm.egGO2ALLEGS)
  genes <- unlist(mget(alleges, org.Mm.egENSEMBL))
  #print(as.vector(genes))
  genes_in_cat <- intersect(as.vector(genes), as.vector(up_degs))
  #print(genes_in_cat)
  gene_sym_in_cat <- as.vector(unlist(mget(unlist(mget(genes_in_cat, org.Mm.egENSEMBL2EG)), org.Mm.egSYMBOL)))
  gene_sym_in_cat_str <- ""
  
  if(length(genes_in_cat) > 0){
    for(j in 1:length(gene_sym_in_cat)){
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j], sep = ',')
    }
  }
  
  go_results$Genes[i] <- gene_sym_in_cat_str
  go_results$no_of_genes[i] <- length(gene_sym_in_cat)
}

write.table(go_results, file = "DESeq2_results/DESeq2result_WT_vs_KO_up_GO.txt", quote=FALSE, sep = '\t', row.names = FALSE)



