# André Vidas Olsen, DCRC, February 2017
rm(list=ls(all=TRUE))


library(biomaRt) # Finding Ensembl IDs
library(edgeR)
library(ggplot2)


# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR MERGING DATA SETS FROM THE LEVEL 1 mRNA PIPELINE
# -----------------------------------------------------------------------------------------------------------------------------
mergeData <- function(path){
  #
  # This funciton merges data sets from the level 1 mRNA pipeline
  #
   
  # find file names
  files <- list.files(path = path,pattern = "\\.tsv$")
  
  #### merge files into data.frame()
  mergedData <- data.frame()
  firstRun <- TRUE
  for(file in files){
    data <- read.table(paste0(path, file), sep = "\t", row.names = 1)
    
    if(firstRun){
      mergedData <- data
      firstRun <- FALSE
    }else{
      mergedData <- cbind(mergedData , data[rownames(mergedData),])
    }
  }
  
  # set colnames to the file names
  colnames(mergedData) <- files
  
  
  return(mergedData)
}



# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR DIFFERENTIAL EXPRESSION ANALYSIS - edgeR,
# ref: M. Robinson et al.: edgeR, differential expression analysis of digital gene expression data. User's Guide.
# -----------------------------------------------------------------------------------------------------------------------------
DEA <- function(data, group, minSamplesInGroups, FC_thresh, FDR_thresh){
  #
  # This function uses glmQLFit() + glmQLTest to process DEA. The functions returns two list with up and down regulated genes respectively.  
  #
  
  y <- DGEList(counts=data,group=group)
  
  
  # filtering: remove lowly expressed genes
  cpmThresh <- 7/min(y$samples$lib.size)*1e6 # set cpmThresh so that it corresponds to 7 (5-10) counts for the minimum library. 
  keep <- rowSums(cpm(y) > cpmThresh) >= minSamplesInGroups # check the lowest library size, 
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  
  # re-calculate the dge object after gene removal
  y <- DGEList(counts=y,group=group)
  
  
  # normalize according to library size
  y <- calcNormFactors(y)
  design <- model.matrix(~group, data = y$samples)
  y <- estimateDisp(y,design)
  
  
  # perform quasi-likelihood F-tests:
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  
  
  # create lists with up/down regulated genes
  tags <- topTags(qlf, n=nrow(y$counts))
  tags <- tags$table
  
  up <- tags[tags$logFC >= FC_thresh & tags$FDR < FDR_thresh, ]
  down <- tags[tags$logFC <= -FC_thresh & tags$FDR < FDR_thresh, ]
  
  DEA_list <- list(up,down)
  names(DEA_list) <- c("up","down")
  
  
  return(DEA_list)
}


# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR REMOVING NON-CODING GENES BASED ON ENSG* ID'S
# -----------------------------------------------------------------------------------------------------------------------------
removeNonCodingGenes <- function(data){
  #
  # This function finds gene names from ENSG* ID's. It then removes all non-coding genes (ENSG* ID's which does not have a gene name)
  #
  
  data_ENSGs <- rownames(data) 
  
  # create Mart object
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  
  # convert ENSG* ID's to gene names and set these in a new column.
  #listOfGenes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=data_ENSGs,mart= mart)
  listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
  
  # merge each gene name to the corresponding ENSG* ID's in the data
  data$geneName <- listOfGenes$hgnc_symbol[match(data_ENSGs, listOfGenes$ensembl_gene_id)]
  
  
  # remove every row without a gene name 
  data <- data[!is.na(data[,"geneName"]),]
  data <- data[ !(data[,"geneName"] == "") ,]
  
  # remove geneName column, since it should not be used in the DEA analysis.
  data <- data[,-ncol(data)]
  
  # remove every mir-gene (start with MIR in the geneName column)
  #data <- data[!grepl("^MIR", data$geneName), ]
  
  
  return(data)
}


# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR FINDING GENE NAMES FROM ENSG* ID'S
# -----------------------------------------------------------------------------------------------------------------------------
findGeneNames <- function(DEA_list){
  #
  # This function finds gene names from ENSG* ID's and merge these names into the input data.frame (DEA_list).
  #
  
  genes_up <- rownames(DEA_list$up)
  genes_down <- rownames(DEA_list$down)
  
  
  # create Mart object
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
   
  # convert ENSG* ID's to gene names and set these in a new column (both formal gene name and full gene name [description]).
  gList_up <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),values=genes_up,mart= mart)
  gList_down <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),values=genes_down,mart= mart) 
  
  
  
  DEA_list$up$geneName <- gList_up$hgnc_symbol[match(rownames(DEA_list$up), gList_up$ensembl_gene_id)]
  DEA_list$up$fullGeneName <- gList_up$description[match(rownames(DEA_list$up), gList_up$ensembl_gene_id)]
  
  DEA_list$down$geneName <- gList_down$hgnc_symbol[match(rownames(DEA_list$down), gList_down$ensembl_gene_id)]
  DEA_list$down$fullGeneName <- gList_down$description[match(rownames(DEA_list$down), gList_down$ensembl_gene_id)]
  
  
  return(DEA_list)
}




# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CREATING PCA OBJECT AND VISUALIZATION
# -----------------------------------------------------------------------------------------------------------------------------
PCAplot <- function(data,group){
  # PCA object
  pca <- prcomp(data, center = TRUE, scale=TRUE)
  rotation <- data.frame(pca$rotation)
  pca$group <- group
  
  # PCA plot
  p <- ggplot(rotation, aes(PC1, PC2)) + geom_text(aes(label= factor(pca$group), colour= factor(pca$group)))
  
  # scree plot
  plot(pca, type = "l")
  #pca$sdev
  
  return(list(pca,p))
}




# -----------------------------------------------------------------------------------------------------------------------------
# LOAD AND CONVERT DATA
# -----------------------------------------------------------------------------------------------------------------------------
# t3_t0 data
path_t3_t0 <- "../../Data/inputData/t3_t0/"
data_t3_t0 <- mergeData(path_t3_t0)

head(data_t3_t0)

# get gene name and remove rows which does not have coding genes
data_t3_t0 <- removeNonCodingGenes(data_t3_t0)

#write.csv(data_t3_t0, file = "../../Data/results/data_t3_t0.csv", quote = FALSE)


# t2_t0 data
path_t2_t0 <- "../../Data/inputData/t2_t0/"
data_t2_t0 <- mergeData(path_t2_t0)

head(data_t2_t0)

# get gene name and remove rows which does not have coding genes
data_t2_t0 <- removeNonCodingGenes(data_t2_t0)

write.csv(data_t2_t0, file = "../../Data/results/data_t2_t0.csv", quote = FALSE)


# t1_t0 data
path_t1_t0 <- "../../Data/inputData/t1_t0/"
data_t1_t0 <- mergeData(path_t1_t0)

head(data_t1_t0)

# get gene name and remove rows which does not have coding genes
data_t1_t0 <- removeNonCodingGenes(data_t1_t0)

write.csv(data_t1_t0, file = "../../Data/results/data_t1_t0.csv", quote = FALSE)



# -----------------------------------------------------------------------------------------------------------------------------
# DEA
# -----------------------------------------------------------------------------------------------------------------------------
# t3_t0 data
colnames(data_t3_t0)
group <- factor(c(1,1,1,2,2,2))
minSamplesInGroups <- 3

DEA_list_t3_t0 <- DEA(data_t3_t0, group, minSamplesInGroups, 1, 0.05)




# t2_t0 data
colnames(data_t2_t0)
group <- factor(c(1,1,1,2,2,2))
minSamplesInGroups <- 3

DEA_list_t2_t0 <- DEA(data_t2_t0, group, minSamplesInGroups, 1, 0.05)




# t1_t0 data
colnames(data_t1_t0)
group <- factor(c(1,1,1,2,2,2))
minSamplesInGroups <- 3

DEA_list_t1_t0 <- DEA(data_t1_t0, group, minSamplesInGroups, 1, 0.05)


# -----------------------------------------------------------------------------------------------------------------------------
# FIND GENE NAMES FROM ENSG ID.
# -----------------------------------------------------------------------------------------------------------------------------
# t3_t0 data
DEA_list_w_geneNames_t3_t0 <- findGeneNames(DEA_list_t3_t0)
DEA_list_w_geneNames_t3_t0

# save up/down data.frames
write.csv(DEA_list_w_geneNames_t3_t0$up, file = "../../Data/results/DE_genes/up_data_t3_t0.csv", quote = TRUE)
write.csv(DEA_list_w_geneNames_t3_t0$down, file = "../../Data/results/DE_genes/down_data_t3_t0.csv", quote = TRUE)




# t2_t0 data
DEA_list_w_geneNames_t2_t0 <- findGeneNames(DEA_list_t2_t0)
DEA_list_w_geneNames_t2_t0

# save up/down data.frames
write.csv(DEA_list_w_geneNames_t2_t0$up, file = "../../Data/results/DE_genes/up_data_t2_t0.csv", quote = TRUE)
write.csv(DEA_list_w_geneNames_t2_t0$down, file = "../../Data/results/DE_genes/down_data_t2_t0.csv", quote = TRUE)




# t1_t0 data
DEA_list_w_geneNames_t1_t0 <- findGeneNames(DEA_list_t1_t0)
DEA_list_w_geneNames_t1_t0

# save up/down data.frames
write.csv(DEA_list_w_geneNames_t1_t0$up, file = "../../Data/results/DE_genes/up_data_t1_t0.csv", quote = TRUE)
write.csv(DEA_list_w_geneNames_t1_t0$down, file = "../../Data/results/DE_genes/down_data_t1_t0.csv", quote = TRUE)



# -----------------------------------------------------------------------------------------------------------------------------
# COMPARISON WITH ELIFE ARTICLE
# -----------------------------------------------------------------------------------------------------------------------------
library(gplots)
# load comparison data (eLife list of genes)
up_genes_elife <- read.table("../../Data/validationData/up_genes_elife.csv", sep = "\t", as.is = TRUE)
down_genes_elife <- read.table("../../Data/validationData/down_genes_elife.csv", sep = "\t", as.is = TRUE)

# t3 vs eLife
up_venn <- venn(list("t3_t0"=DEA_list_w_geneNames_t3_t0$up$geneName, "eLife"=up_genes_elife$V1))
title("Venn diagram between t3_t0 and eLife \n (up-regulated genes)")

intersect_up_t3_eLife <- intersect(DEA_list_w_geneNames_t3_t0$up$geneName, up_genes_elife$V1)
write.csv(intersect_up_t3_eLife, file = "../../Data/results/overlap/intersect_up_t3_eLife.csv", quote = FALSE)


down_venn <- venn(list("t3_t0"=DEA_list_w_geneNames_t3_t0$down$geneName, "eLife"=down_genes_elife$V1))
title("Venn diagram between t3_t0 and eLife \n (down-regulated genes)")

intersect_down_t3_eLife <- intersect(DEA_list_w_geneNames_t3_t0$down$geneName, down_genes_elife$V1)
write.csv(intersect_down_t3_eLife, file = "../../Data/results/overlap/intersect_down_t3_eLife.csv", quote = FALSE)


# t2 vs eLife
up_venn <- venn(list("t2_t0"=DEA_list_w_geneNames_t2_t0$up$geneName, "eLife"=up_genes_elife$V1))
title("Venn diagram between t2_t0 and eLife \n (up-regulated genes)")

intersect_up_t2_eLife <- intersect(DEA_list_w_geneNames_t2_t0$up$geneName, up_genes_elife$V1)
write.csv(intersect_up_t2_eLife, file = "../../Data/results/overlap/intersect_up_t2_eLife.csv", quote = FALSE)


down_venn <- venn(list("t2_t0"=DEA_list_w_geneNames_t2_t0$down$geneName, "eLife"=down_genes_elife$V1))
title("Venn diagram between t2_t0 and eLife \n (down-regulated genes)")

intersect_down_t2_eLife <- intersect(DEA_list_w_geneNames_t2_t0$down$geneName, down_genes_elife$V1)
write.csv(intersect_down_t2_eLife, file = "../../Data/results/overlap/intersect_down_t2_eLife.csv", quote = FALSE)


# t1 vs eLife
up_venn <- venn(list("t1_t0"=DEA_list_w_geneNames_t1_t0$up$geneName, "eLife"=up_genes_elife$V1))
title("Venn diagram between t1_t0 and eLife \n (up-regulated genes)")

intersect_up_t1_eLife <- intersect(DEA_list_w_geneNames_t1_t0$up$geneName, up_genes_elife$V1)
write.csv(intersect_up_t1_eLife, file = "../../Data/results/overlap/intersect_up_t1_eLife.csv", quote = FALSE)


down_venn <- venn(list("t1_t0"=DEA_list_w_geneNames_t1_t0$down$geneName, "eLife"=down_genes_elife$V1))
title("Venn diagram between t1_t0 and eLife \n (down-regulated genes)")

intersect_down_t1_eLife <- intersect(DEA_list_w_geneNames_t1_t0$down$geneName, down_genes_elife$V1)
write.csv(intersect_down_t1_eLife, file = "../../Data/results/overlap/intersect_down_t1_eLife.csv", quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------
# COMPARISON BETWEEN TIME POINTS
# -----------------------------------------------------------------------------------------------------------------------------



# t2 vs t1
up_venn <- venn(list("t2_t0"=DEA_list_w_geneNames_t2_t0$up$geneName, "t1_t0"=DEA_list_w_geneNames_t1_t0$up$geneName))
title("Venn diagram between t2_t0 and t1_t0 \n (up-regulated genes)")

intersect_up_t2_t1 <- intersect(DEA_list_w_geneNames_t2_t0$up$geneName, DEA_list_w_geneNames_t1_t0$up$geneName)
write.csv(intersect_up_t2_t1, file = "../../Data/results/overlap/intersect_up_t2_t1.csv", quote = FALSE)


down_venn <- venn(list("t2_t0"=DEA_list_w_geneNames_t2_t0$down$geneName, "t1_t0"=DEA_list_w_geneNames_t1_t0$down$geneName))
title("Venn diagram between t2_t0 and t1_t0 \n (down-regulated genes)")

intersect_down_t2_t1 <- intersect(DEA_list_w_geneNames_t2_t0$down$geneName, DEA_list_w_geneNames_t1_t0$down$geneName)
write.csv(intersect_down_t2_t1, file = "../../Data/results/overlap/intersect_down_t2_t1.csv", quote = FALSE)


# t3 vs t2
up_venn <- venn(list("t3_t0"=DEA_list_w_geneNames_t3_t0$up$geneName, "t2_t0"=DEA_list_w_geneNames_t2_t0$up$geneName))
title("Venn diagram between t3_t0 and t2_t0 \n (up-regulated genes)")

intersect_up_t3_t2 <- intersect(DEA_list_w_geneNames_t3_t0$up$geneName, DEA_list_w_geneNames_t2_t0$up$geneName)
write.csv(intersect_up_t3_t2, file = "../../Data/results/overlap/intersect_up_t3_t2.csv", quote = FALSE)


down_venn <- venn(list("t3_t0"=DEA_list_w_geneNames_t3_t0$down$geneName, "t2_t0"=DEA_list_w_geneNames_t2_t0$down$geneName))
title("Venn diagram between t3_t0 and t2_t0 \n (down-regulated genes)")

intersect_down_t3_t2 <- intersect(DEA_list_w_geneNames_t3_t0$down$geneName, DEA_list_w_geneNames_t2_t0$down$geneName)
write.csv(intersect_down_t3_t2, file = "../../Data/results/overlap/intersect_down_t3_t2.csv", quote = FALSE)


# t3 vs t1
up_venn <- venn(list("t3_t0"=DEA_list_w_geneNames_t3_t0$up$geneName, "t1_t0"=DEA_list_w_geneNames_t1_t0$up$geneName))
title("Venn diagram between t3_t0 and t1_t0 \n (up-regulated genes)")

intersect_up_t3_t1 <- intersect(DEA_list_w_geneNames_t3_t0$up$geneName, DEA_list_w_geneNames_t1_t0$up$geneName)
write.csv(intersect_up_t3_t1, file = "../../Data/results/overlap/intersect_up_t3_t1.csv", quote = FALSE)


down_venn <- venn(list("t3_t0"=DEA_list_w_geneNames_t3_t0$down$geneName, "t1_t0"=DEA_list_w_geneNames_t1_t0$down$geneName))
title("Venn diagram between t3_t0 and t1_t0 \n (down-regulated genes)")

intersect_down_t3_t1 <- intersect(DEA_list_w_geneNames_t3_t0$down$geneName, DEA_list_w_geneNames_t1_t0$down$geneName)
write.csv(intersect_down_t3_t1, file = "../../Data/results/overlap/intersect_down_t3_t1.csv", quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------
# VISUALIZATION - VENN DIAGRAM FOR OVERLAPPING RESULTS BETWEEN DATA AND RAW DATA 
# -----------------------------------------------------------------------------------------------------------------------------
library(gplots)


# up-regulated genes
#up_venn <- venn(list("Data"=DEA_list_w_geneNames$up$geneName, "Raw Data"=DEA_list_w_geneNames_raw$up$geneName)) + title("Venn diagram between data and raw data \n (up-regulated genes)")
      

# down-regulated genes
#down_venn <- venn(list("Data"=DEA_list_w_geneNames$down$geneName, "Raw Data"=DEA_list_w_geneNames_raw$down$geneName)) + title("Venn diagram between data and raw data \n (down-regulated genes)")


# -----------------------------------------------------------------------------------------------------------------------------
# VISUALIZATION - PCA-PLOT 
# -----------------------------------------------------------------------------------------------------------------------------





# t3_t0 data before the DEA
head(data_t3_t0)

PCA <- PCAplot(data_t3_t0,c(0,0,0,3,3,3))
PCA[2]







# t3_t0 data containing only the differentially expressed genes
data_t3_t0_DE_genes <- data_t3_t0[match(c(rownames(DEA_list_w_geneNames_t3_t0$up),rownames(DEA_list_w_geneNames_t3_t0$down)), rownames(data_t3_t0)),]
head(data_t3_t0_DE_genes)

PCA <- PCAplot(data_t3_t0_DE_genes,c(0,0,0,3,3,3))
PCA[2]


