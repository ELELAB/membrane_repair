rm(list=ls())
library(pathview)
library(biomaRt)
library(DOSE)
library(clusterProfiler)

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
  entrezID <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'hgnc_symbol', values = data$gene , mart = mart)
  
  
  # insert fc into entrezID data set, so only genes with entrezID's exits in the data set.
  entrezID$foldChange <- data$logFC[match(entrezID$hgnc_symbol, data$gene)]
  
  
  # remove NA's from entrezID data set
  entrezID <- entrezID[!is.na(entrezID$entrezgene),]
  
  #print(entrezID)
  
  foldChange <- entrezID$foldChange 
  
  # insert entrezID as rownames
  names(foldChange) <- entrezID$entrezgene
  
  return(foldChange)
}



t3_t0_up <- read.csv(file = "../../Data/results/DE_genes/up_data_t3_t0.csv", header = TRUE, row.names = 1)
t3_t0_down <- read.csv(file = "../../Data/results/DE_genes/down_data_t3_t0.csv", header = TRUE, row.names = 1)
t3_t0 <- rbind(t3_t0_up,t3_t0_down)
head(t3_t0)

t2_t0_up <- read.csv(file = "../../Data/results/DE_genes/up_data_t2_t0.csv", header = TRUE, row.names = 1)
t2_t0_down <- read.csv(file = "../../Data/results/DE_genes/down_data_t2_t0.csv", header = TRUE, row.names = 1)
t2_t0 <- rbind(t2_t0_up,t2_t0_down)
head(t2_t0)

t1_t0_up <- read.csv(file = "../../Data/results/DE_genes/up_data_t1_t0.csv", header = TRUE, row.names = 1)
t1_t0_down <- read.csv(file = "../../Data/results/DE_genes/down_data_t1_t0.csv", header = TRUE, row.names = 1)
t1_t0 <- rbind(t1_t0_up,t1_t0_down)
head(t1_t0)


# find list of ordered enriched pathways
de_t3_t0 <- names(convertToPathviewFormat(t3_t0))
de_t2_t0 <- names(convertToPathviewFormat(t2_t0))
de_t1_t0 <- names(convertToPathviewFormat(t1_t0))

kk_t3_t0 <- enrichKEGG(de_t3_t0, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
kk_t3_t0_list <- summary(kk_t3_t0)$ID
write.csv(summary(kk_t3_t0), file = "pathwayList_t3_t0.csv", quote = FALSE)

kk_t2_t0 <- enrichKEGG(de_t2_t0, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
kk_t2_t0_list <- summary(kk_t2_t0)$ID
write.csv(summary(kk_t2_t0), file = "pathwayList_t2_t0.csv", quote = FALSE)

kk_t1_t0 <- enrichKEGG(de_t1_t0, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
kk_t1_t0_list <- summary(kk_t1_t0)$ID
write.csv(summary(kk_t1_t0), file = "pathwayList_t1_t0.csv", quote = FALSE)


# enrichment analysis

# t3_t0
foldChange <- convertToPathviewFormat(t3_t0)
for(element in kk_t3_t0_list){
  print(unlist(strsplit(element,split = "hsa"))[2])
  pv.out <- pathview(gene.data = foldChange, pathway.id = unlist(strsplit(element,split = "hsa"))[2], species = "hsa", out.suffix = "t3_t0", kegg.native = T)
}



# t2_t0
foldChange <- convertToPathviewFormat(t2_t0)
for(element in kk_t2_t0_list){
  print(unlist(strsplit(element,split = "hsa"))[2])
  pv.out <- pathview(gene.data = foldChange, pathway.id = unlist(strsplit(element,split = "hsa"))[2], species = "hsa", out.suffix = "t2_t0", kegg.native = T)
}


# t1_t0
foldChange <- convertToPathviewFormat(t1_t0)
for(element in kk_t1_t0_list){
  print(unlist(strsplit(element,split = "hsa"))[2])
  pv.out <- pathview(gene.data = foldChange, pathway.id = unlist(strsplit(element,split = "hsa"))[2], species = "hsa", out.suffix = "t1_t0", kegg.native = T)
}




