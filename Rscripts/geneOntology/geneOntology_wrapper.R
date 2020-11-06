rm(list=ls(all=TRUE))
source("geneOntology_functions.R")



# ---------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------
findGeneList <- function(data){
  #
  # This function finds the gene lists of all rows in a given count data set (rows: protein coding ENSG*, cols: samples)
  # Requires the following packages: biomaRt
  library(biomaRt) # Finding Ensembl IDs
  
  data_ENSGs <- rownames(data) 
  
  # create Mart object
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
  
  # merge each gene name to the corresponding ENSG* ID's in the data
  listOfGenes <- listOfGenes$hgnc_symbol[match(data_ENSGs, listOfGenes$ensembl_gene_id)]
  
  return(listOfGenes)
}

# ---------------------------------------------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------------------------------------------
#### load HUGO
load("HUGO.RData")

#### gene universe

# t3_t0
t3_t0 <- read.csv(file = "../../Data/results/data_t3_t0.csv", row.names = 1)
geneUni_t3_t0 <- findGeneList(t3_t0)


# t2_t0
t2_t0 <- read.csv(file = "../../Data/results/data_t2_t0.csv", row.names = 1)
geneUni_t2_t0 <- findGeneList(t2_t0)


# t1_t0
t1_t0 <- read.csv(file = "../../Data/results/data_t1_t0.csv", row.names = 1)
geneUni_t1_t0 <- findGeneList(t1_t0)



#### DE gene lists

# t3_t0
up_t3_t0 <- read.csv(file = "../../Data/results/DE_genes/up_data_t3_t0.csv", row.names = 1, as.is = TRUE)
up_t3_t0_geneList <- up_t3_t0$geneName


down_t3_t0 <- read.csv(file = "../../Data/results/DE_genes/down_data_t3_t0.csv", row.names = 1, as.is = TRUE)
down_t3_t0_geneList <- down_t3_t0$geneName

# t2_t0
up_t2_t0 <- read.csv(file = "../../Data/results/DE_genes/up_data_t2_t0.csv", row.names = 1, as.is = TRUE)
up_t2_t0_geneList <- up_t2_t0$geneName

down_t2_t0 <- read.csv(file = "../../Data/results/DE_genes/down_data_t2_t0.csv", row.names = 1, as.is = TRUE)
down_t2_t0_geneList <- down_t2_t0$geneName

# t1_t0
up_t1_t0 <- read.csv(file = "../../Data/results/DE_genes/up_data_t1_t0.csv", row.names = 1, as.is = TRUE)
up_t1_t0_geneList <- up_t1_t0$geneName

down_t1_t0 <- read.csv(file = "../../Data/results/DE_genes/down_data_t1_t0.csv", row.names = 1, as.is = TRUE)
down_t1_t0_geneList <- down_t1_t0$geneName


# ---------------------------------------------------------------------------------------------------------------
# Find GO terms for BP
# ---------------------------------------------------------------------------------------------------------------
t3_t0_GO_up <- TOPGO("BP", geneUni_t3_t0, up_t3_t0_geneList, HUGO, 30, "t3_t0_GO_up")
write.csv(t3_t0_GO_up, file = "../../Data/results/GO/BP/t3_t0_GO_up.csv", quote = TRUE)



t3_t0_GO_down <- TOPGO("BP", geneUni_t3_t0, down_t3_t0_geneList, HUGO, 30, "t3_t0_GO_down")
write.csv(t3_t0_GO_down, file = "../../Data/results/GO/BP/t3_t0_GO_down.csv", quote = TRUE)



t2_t0_GO_up <- TOPGO("BP", geneUni_t2_t0, up_t2_t0_geneList, HUGO, 30, "t2_t0_GO_up")
write.csv(t2_t0_GO_up, file = "../../Data/results/GO/BP/t2_t0_GO_up.csv", quote = TRUE)

t2_t0_GO_down <- TOPGO("BP", geneUni_t2_t0, down_t2_t0_geneList, HUGO, 30, "t2_t0_GO_down")
write.csv(t2_t0_GO_down, file = "../../Data/results/GO/BP/t2_t0_GO_down.csv", quote = TRUE)



t1_t0_GO_up <- TOPGO("BP", geneUni_t1_t0, up_t1_t0_geneList, HUGO, 30, "t1_t0_GO_up")
write.csv(t1_t0_GO_up, file = "../../Data/results/GO/BP/t1_t0_GO_up.csv", quote = TRUE)

t1_t0_GO_down <- TOPGO("BP", geneUni_t1_t0, down_t1_t0_geneList, HUGO, 30, "t1_t0_GO_down")
write.csv(t1_t0_GO_down, file = "../../Data/results/GO/BP/t1_t0_GO_down.csv", quote = TRUE)



# ---------------------------------------------------------------------------------------------------------------
# Find GO terms for CC
# ---------------------------------------------------------------------------------------------------------------
t3_t0_GO_up <- TOPGO("CC", geneUni_t3_t0, up_t3_t0_geneList, HUGO, 30, "t3_t0_GO_up")
write.csv(t3_t0_GO_up, file = "../../Data/results/GO/CC/t3_t0_GO_up.csv", quote = TRUE)

t3_t0_GO_down <- TOPGO("CC", geneUni_t3_t0, down_t3_t0_geneList, HUGO, 30, "t3_t0_GO_down")
write.csv(t3_t0_GO_down, file = "../../Data/results/GO/CC/t3_t0_GO_down.csv", quote = TRUE)



t2_t0_GO_up <- TOPGO("CC", geneUni_t2_t0, up_t2_t0_geneList, HUGO, 30, "t2_t0_GO_up")
write.csv(t2_t0_GO_up, file = "../../Data/results/GO/CC/t2_t0_GO_up.csv", quote = TRUE)

t2_t0_GO_down <- TOPGO("CC", geneUni_t2_t0, down_t2_t0_geneList, HUGO, 30, "t2_t0_GO_down")
write.csv(t2_t0_GO_down, file = "../../Data/results/GO/CC/t2_t0_GO_down.csv", quote = TRUE)



t1_t0_GO_up <- TOPGO("CC", geneUni_t1_t0, up_t1_t0_geneList, HUGO, 30, "t1_t0_GO_up")
write.csv(t1_t0_GO_up, file = "../../Data/results/GO/CC/t1_t0_GO_up.csv", quote = TRUE)

t1_t0_GO_down <- TOPGO("CC", geneUni_t1_t0, down_t1_t0_geneList, HUGO, 30, "t1_t0_GO_down")
write.csv(t1_t0_GO_down, file = "../../Data/results/GO/CC/t1_t0_GO_down.csv", quote = TRUE)


# ---------------------------------------------------------------------------------------------------------------
# Find GO terms for MF
# ---------------------------------------------------------------------------------------------------------------
t3_t0_GO_up <- TOPGO("MF", geneUni_t3_t0, up_t3_t0_geneList, HUGO, 30, "t3_t0_GO_up")
write.csv(t3_t0_GO_up, file = "../../Data/results/GO/MF/t3_t0_GO_up.csv", quote = TRUE)

t3_t0_GO_down <- TOPGO("MF", geneUni_t3_t0, down_t3_t0_geneList, HUGO, 30, "t3_t0_GO_down")
write.csv(t3_t0_GO_down, file = "../../Data/results/GO/MF/t3_t0_GO_down.csv", quote = TRUE)



t2_t0_GO_up <- TOPGO("MF", geneUni_t2_t0, up_t2_t0_geneList, HUGO, 30, "t2_t0_GO_up")
write.csv(t2_t0_GO_up, file = "../../Data/results/GO/MF/t2_t0_GO_up.csv", quote = TRUE)

t2_t0_GO_down <- TOPGO("MF", geneUni_t2_t0, down_t2_t0_geneList, HUGO, 30, "t2_t0_GO_down")
write.csv(t2_t0_GO_down, file = "../../Data/results/GO/MF/t2_t0_GO_down.csv", quote = TRUE)



t1_t0_GO_up <- TOPGO("MF", geneUni_t1_t0, up_t1_t0_geneList, HUGO, 30, "t1_t0_GO_up")
write.csv(t1_t0_GO_up, file = "../../Data/results/GO/MF/t1_t0_GO_up.csv", quote = TRUE)

t1_t0_GO_down <- TOPGO("MF", geneUni_t1_t0, down_t1_t0_geneList, HUGO, 30, "t1_t0_GO_down")
write.csv(t1_t0_GO_down, file = "../../Data/results/GO/MF/t1_t0_GO_down.csv", quote = TRUE)





