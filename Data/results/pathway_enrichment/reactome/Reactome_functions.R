
library(ReactomePA)
library(biomaRt)

# --------------------------------------------------------------------------------------
# PATHWAY ENRICHMENT ANALYSIS

# background: dataframe that includes the universe gene entrez-IDs and the corrisponding
# gene symbols
# my.genes: genes selected as DE
# my.name: the name of the output file
# my.plot: logical (TRUE if you want an output plot)
# value: p-value cutoff (default: 0.05)
# --------------------------------------------------------------------------------------


enrich_pathway <- function(background, my.genes, my.name, my.plot, value=0.05) {
  entrez.data <- background[background$geneName %in% my.genes, ]
  my.pathway <- enrichPathway(as.character(entrez.data$entrez), organism = "human", pvalueCutoff = value, pAdjustMethod = "fdr", qvalueCutoff = value, as.character(background$entrez), minGSSize = 3, readable = T)
  if (my.plot == TRUE) {
    png(filename = paste0(my.name,".png"), height = 1200, width = 1200)
    cnetplot(my.pathway, categorySize="pvalue")
    dev.off()
  }
  return(my.pathway)
}

#------------------------------------------------------------
# convert the total gene-names to entrezIDs

# universe_genes: the expression matrix with all genes
#----------------------------------------------------------

conversion.background <- function(universe_genes){
  total_genes <- rownames(universe_genes)
  #covert the gene names (universal genes) into entrez IDs
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  background <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'hgnc_symbol',
                      values = total_genes , mart = ensembl)
  colnames(background) <- c("geneName","entrez")
  return (background)
}


convertToPathviewFormat <- function(data){
  #
  # This function convert data consisting of gene name (gene col.) and fold change values (foldChangeUp/foldChangeDown) to the data format used by pathview
  #
  
  # merge dataUp and dataDown together
  #dataDown$foldChangeDown <- dataDown$foldChangeDown*(-1) # convert FC down to negative values, so they fit with FC up data.
  #data <- data.frame(gene = c(dataUp$gene[1:10], dataDown$gene[1:10]), foldChange = c(dataUp$foldChangeUp[1:10], dataDown$foldChangeDown[1:10]))
  
  
  
  # convert to entrez IDs
  mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  # to get the corresponding entrez IDs of my DE genes symbols
  entrezID <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'ensembl_gene_id', values = rownames(data) , mart = mart)
  
  # remove NA's from entrezID data set
  entrezID <- entrezID[!is.na(entrezID$entrezgene),]
  
  # change colnames
  colnames(entrezID) <- c("geneName","entrez")
  
  return(entrezID)
}


