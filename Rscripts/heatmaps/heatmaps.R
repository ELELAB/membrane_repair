rm(list=ls(all=TRUE))
library(biomaRt) # Finding Ensembl IDs
library(edgeR)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(latex2exp)


#### get data. Controls: 220, 221, 222
# expression data
#t1_t0 <- read.csv("../../Data/results/data_t1_t0.csv", row.names = 1 , stringsAsFactors = FALSE)
#t2_t0 <- read.csv("../../Data/results/data_t2_t0.csv", row.names = 1, stringsAsFactors = FALSE)
#t3_t0 <- read.csv("../../Data/results/data_t3_t0.csv", row.names = 1, stringsAsFactors = FALSE)


# DE genes
DEA_t1_down <- read.csv("../../Data/results/DE_genes/down_data_t1_t0.csv", stringsAsFactors = FALSE)
DEA_t1_up <- read.csv("../../Data/results/DE_genes/up_data_t1_t0.csv", stringsAsFactors = FALSE)
DEA_t1 <- rbind(DEA_t1_down, DEA_t1_up)  
DEA_t2_down <- read.csv("../../Data/results/DE_genes/down_data_t2_t0.csv", stringsAsFactors = FALSE)
DEA_t2_up <- read.csv("../../Data/results/DE_genes/up_data_t2_t0.csv", stringsAsFactors = FALSE)
DEA_t2 <- rbind(DEA_t2_down, DEA_t2_up)
DEA_t3_down <- read.csv("../../Data/results/DE_genes/down_data_t3_t0.csv", stringsAsFactors = FALSE)
DEA_t3_up <- read.csv("../../Data/results/DE_genes/up_data_t3_t0.csv", stringsAsFactors = FALSE)
DEA_t3 <- rbind(DEA_t3_down, DEA_t3_up)


# test purposes
#DEA_t2_down[, c("geneName", "logFC")]
#DEA_t1[DEA_t1$geneName == "C1orf111", ]
#DEA_t2[DEA_t2$geneName == "C1orf111", ]
#dim(DEA_t3)
#DEA_t3 <- DEA_t3[!DEA_t3$geneName == "C1orf111", ]
#dim(DEA_t3)


# DE and non genes from DEA (p = 1, FC = 0)
DEA_t1_down_all <- read.csv("../../Data/results/DE_genes/heatmaps_down_data_t1_t0.csv", stringsAsFactors = FALSE)
DEA_t1_up_all <- read.csv("../../Data/results/DE_genes/heatmaps_up_data_t1_t0.csv", stringsAsFactors = FALSE)
DEA_t1_all <- rbind(DEA_t1_down_all, DEA_t1_up_all)  
DEA_t2_down_all <- read.csv("../../Data/results/DE_genes/heatmaps_down_data_t2_t0.csv", stringsAsFactors = FALSE)
DEA_t2_up_all <- read.csv("../../Data/results/DE_genes/heatmaps_up_data_t2_t0.csv", stringsAsFactors = FALSE)
DEA_t2_all <- rbind(DEA_t2_down_all, DEA_t2_up_all)
DEA_t3_down_all <- read.csv("../../Data/results/DE_genes/heatmaps_down_data_t3_t0.csv", stringsAsFactors = FALSE)
DEA_t3_up_all <- read.csv("../../Data/results/DE_genes/heatmaps_up_data_t3_t0.csv", stringsAsFactors = FALSE)
DEA_t3_all <- rbind(DEA_t3_down_all, DEA_t3_up_all)

DEA_t1_all[DEA_t1_all$geneName == "C1orf111", ]
DEA_t2_all[DEA_t2_all$geneName == "C1orf111", ]
DEA_t3_all[DEA_t3_all$geneName == "C1orf111", ]

# find pool of all regulated genes in DEA_t1-t3
DE_gene_pool <- unique(c(DEA_t1$geneName, DEA_t2$geneName, DEA_t3$geneName))



# extract rows from DEA_*_all DE and non-DE gene results to obtain FC values for all genes (including non-DE genes)
DEA_t1_all <- DEA_t1_all[DEA_t1_all$geneName %in% DE_gene_pool, ]
DEA_t2_all <- DEA_t2_all[DEA_t2_all$geneName %in% DE_gene_pool, ]
DEA_t3_all <- DEA_t3_all[DEA_t3_all$geneName %in% DE_gene_pool, ]


# make heatmap matrix
data <- matrix(nrow = length(DE_gene_pool), ncol = 3, dimnames = list(DE_gene_pool, c("t1","t2","t3")))


data[,"t1"] <- DEA_t1_all$logFC[match(rownames(data), DEA_t1_all$geneName)]
data[,"t2"] <- DEA_t2_all$logFC[match(rownames(data), DEA_t2_all$geneName)]
data[,"t3"] <- DEA_t3_all$logFC[match(rownames(data), DEA_t3_all$geneName)]



#### extract only DE genes from t1_t0 - t3_ t0 expression data
# create Mart object
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# convert ENSG* ID's to gene names and set these in a new column (formal gene name).
#t1_genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(t1_t0),mart= mart)
#t2_genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(t2_t0),mart= mart)
#t3_genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(t3_t0),mart= mart)
#t1_t3_genes <- rbind(t1_genes, t2_genes, t3_genes)

#DE_genes_ensg_t1 <- unique(t1_genes[t1_genes$hgnc_symbol %in% DE_gene_pool, "ensembl_gene_id"])
#DE_genes_ensg_t2 <- unique(t2_genes[t2_genes$hgnc_symbol %in% DE_gene_pool, "ensembl_gene_id"])
#DE_genes_ensg_t3 <- unique(t3_genes[t3_genes$hgnc_symbol %in% DE_gene_pool, "ensembl_gene_id"])


# extract expression data to only include DE genes
#t1_t0_only_DE_genes <- t1_t0[rownames(t1_t0) %in% DE_genes_ensg_t1,]
#t2_t0_only_DE_genes <- t2_t0[rownames(t2_t0) %in% DE_genes_ensg_t2,]
#t3_t0_only_DE_genes <- t3_t0[rownames(t3_t0) %in% DE_genes_ensg_t3,]



#colnames(t1_t0_only_DE_genes)
#colnames(t2_t0_only_DE_genes)
#colnames(t3_t0_only_DE_genes)

#### transform data values to logFC
#normal <- data[1] # background values to be used for calucation of logFC
#data <- data[,2:ncol(data)] # remove the normal (background) column from the data set
#dim(data)


# calculate logFC for all columns in data
#for(i in 1:ncol(data)){
#  data[,i] <- log2(data[,i]/normal)
#}
#head(data)





#### samples - samples
dataCor <- cor(data, method = "spearman")
head(dataCor)
dataCorMatrix <- data.matrix(dataCor)
Heatmap(dataCorMatrix, name = "Spearman cor.", column_title = "Samples", column_title_side = "bottom", row_title = "Samples", row_title_side = "right")


# classes - classes
dataTransposed <- t(data)
dataCor <- cor(dataTransposed, method = "spearman")
head(dataCor)
dataCorMatrix <- data.matrix(dataCor)
Heatmap(dataCorMatrix, name = "Spearman cor.", column_title = "Classes", column_title_side = "bottom", row_title = "Classes", row_title_side = "right")




# classes - samples

# change -Inf to something very small
dataMatrix <- data.matrix(data)
head(dataMatrix)
dataMatrix <- t(dataMatrix)
dataMatrix[is.na(dataMatrix)] <- 0
head(dataMatrix)

min(dataMatrix)
max(dataMatrix)

bmp(filename="fullHeatmap.bmp", width = 3200, height = 600)
Heatmap(dataMatrix, name = "LogFC", column_title = "Samples", column_title_side = "bottom", row_title = "Genes", row_title_side = "left", row_names_side = "left", na_col = "black", col = colorRamp2(c(min(min(dataMatrix), -max(dataMatrix)), 0, max(-min(dataMatrix), max(dataMatrix))), c("blue", "white", "red")), width = unit(100, "cm"), row_dend_side = "right") #, row_names_gp = gpar(fontsize = 0.1))#, show_row_names = FALSE)
dev.off()
?Heatmap()









nba.m <- melt(dataMatrix)
head(nba.m)
(p <- ggplot(nba.m, aes(X1, X2)) + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue"))




nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
nba$Name <- with(nba, reorder(Name, PTS))
nba[nba$MIN > 38, "MIN"] <- NA
head(nba)
library(ggplot2)
library(reshape)
library(plyr)
library(scales)
nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform, rescale = rescale(value))
(p <- ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue"))
